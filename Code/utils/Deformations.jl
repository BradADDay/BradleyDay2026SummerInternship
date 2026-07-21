using CSV, DataFrames, Interpolations, Gradus

const DeformationBoundsTable = DataFrame(CSV.File("Code/utils/DeformationBounds.csv"))

# ====================================================================================
# Lines
# ====================================================================================

# Get the gradient of a line between two coordinates
function Gradient(coord1, coord2)
    diff = coord1 .- coord2
    diff[2]/diff[1]
end

# Get the intercept of a line with a given gradient that passes through a given point
Intercept(grad, coord) = coord[2] - grad * coord[1]

# Line struct
struct BoundLine{T<:Real}
    "Gradient"
    m::T
    "Intercept"
    c::T
end

# Constructor from two coordinates
function BoundLine(coord1::AbstractArray, coord2::AbstractArray)
    m = Gradient(coord1, coord2)
    c = Intercept(m, coord1)
    BoundLine(m, c)
end

# Get a point on a line
LineEquation(x, line::BoundLine) = line.m * x .+ line.c

# Find the intersection of two lines
Intersection(line1::BoundLine, line2::BoundLine) = (line2.c - line1.c) / (line1.m - line2.m)

# Get the boundary line parameters for a given black hole spin from the bounds table
function GetLineParams(a, df=DeformationBoundsTable)

    # Read the file and separate the spins
    as = df.a
    data = select!(copy(df), Not(:a))

    len = length(names(data))

    # Empty array for storage
    out = Array{Float64}(undef, len)

    # Interpolating the boundary table to get the line parameters
    for (i,name) in enumerate(String.(names(data)))
        interp = linear_interpolation(as, data[!, name])
        out[i] = interp(a)
    end

    out, len

end

# Get the boundary lines for a given black hole spin
function GetLines(a, df=DeformationBoundsTable)
    """
    Read the boundary line definitions, interpolate and return line objects for each of them
    """
    
    out, len = GetLineParams(a, df)

    returns = Array{BoundLine}(undef, Int(round(len/2)))

    # Creating the line objects
    for i in 1:Int(round(len/2))
        returns[i] = BoundLine(out[2*i-1], out[2*i])
    end

    returns
end

# ====================================================================================
# Validity
# ====================================================================================

function ValidityCheckISCO(m)
    """
    Checking if the metric is valid
        returns the ISCO if so
        returns an error value if not
    """
    if is_naked_singularity(m)
        # Naked singularity
        return -2
    elseif ((m.α13 < Constraints(m.a)) | (m.ϵ3 < Constraints(m.a)))
        # Abnormal exterior region
        return -3
    elseif is_no_isco(m)
        # No ISCO
        return -1
    else
        # No Abnormalities
        return Gradus.isco(m)
    end

end

function IsValid(ϵ3, α13, a)
    """
    Check if a parameter combination is valid by explicitly checking the resulting metric
    """
    m = JohannsenMetric(a=a, ϵ3=ϵ3, α13=α13)
    
    if is_naked_singularity(m) | is_no_isco(m) | (ϵ3 < Constraints(a)) | (α13 < Constraints(a)) | (abs(a) > 0.998)
        return false
    else 
        return true
    end

end

function QuickIsValid(ϵ3, α13, a, df=DeformationBoundsTable; sgns=[1,1,-1,-1,-1,-1])
    """
    Checking if a parameter combination is valid by checking if it falls within the pre-defined bounds
    """

    if (ϵ3 < Constraints(a)) | (α13 < Constraints(a)) | (abs(a) > 0.998)
        return false
    end
    
    conds = GetValidityConditions(ϵ3, α13, a, df; sgns=[1,1,-1,-1,-1,-1])
    
    # Checking if the parameter combination is within any of the 4 forbidden regions
    if conds[1] | conds[2] | conds[3] | conds[4]
        return false
    else
        return true
    end

end

function GetValidityConditions(ϵ3, α13, a, df=DeformationBoundsTable; sgns=[1,1,-1,-1,-1,-1])

    # Getting the boundary lines and defining conditions
    lines = GetLines(a, df)
    conds = fill(true, length(lines))

    # Checking if the parameter combination is within the forbidden region for every line
    for i in eachindex(lines)
        conds[i] = sgns[i]*α13 > sgns[i]*(LineEquation(ϵ3, lines[i]))
    end

    conds = ((conds[1] & conds[2]), (conds[3] & conds[4]), conds[5], conds[6])

    return conds

end

function is_no_isco(m)
    """
    Check if there is a valid ISCO for the input metric
    """
    try 
        Gradus.isco(m)
        return false
    catch
        return true
    end
end

# Defining a dispatch for a set of parameters as input
is_no_isco(ϵ3, α13, a) = is_no_isco(JohannsenMetric(ϵ3=ϵ3, α13=α13, a=a))

# The lower limit for the deformation parameters α13, ϵ3 for a given spin
Constraints(a) = -(1+sqrt(1 - a^2))^3

function IsValidFit(fit)
    # Check if a SpectralFitting FitResult is valid
    u=fit.u
    IsValid(u[7], u[6], u[3])
end

# ====================================================================================
# Corrections
# ====================================================================================

function FindNearestSafePoint(point, a; maxIter=1000)
    
    combo = "$a, $point"

    # Moving the coordinate to be within valid parameter space
    point .= ReturnToConstraints.(point, [a]) .+ 0.1

    # Getting the boundary lines
    lines = GetLines(a)

    fixed=false

    iter = 0

    conds=nothing

    for i in 1:maxIter

        # A catch for if the point gets stuck between the top and bottom triangle
        if iter > 3
            point = [point[1] + 0.2, point[2]]
        end

        if (any(point .< Constraints(a)))
            point = ReturnToConstraints.(point, [a])
        end

        # Finding which of the 4 regions the point is in
        conds = GetValidityConditions(point..., a)

        if conds[1] & conds[2]
            # If the point is in both the top and bottom triangle, 
            # moving it to the top one and finding the nearest line
            point = [point[1], point[2] + 1]
            line=NearestLine(point, lines[1:4], a)
        elseif conds[1] | conds[2]
            # If the point is in either the top or bottom triangle,
            # finding the nearest line
            line=NearestLine(point, lines[1:4], a)
            iter +=1
        elseif conds[3]
            # Bottom wedge
            line = lines[5]
        elseif conds[4]
            # Unpredictable region in lower left
            line = lines[6]
        elseif point[1] > 10
            # Checking the point doesn't exceed the set region
            point = [point[1]-0.5, point[2]]
            line = lines[2]
        else
            # Breaking the loop if the point is now in a safe region
            fixed=true
            break
        end

        # Moving the point to the closest point on the closest boundary it is breaching
        point = GetNewPoint(point, line)
 
    end

    if fixed
        return round.(point; digits=2)
    else
        throw(error("Combination $combo could not be corrected! Ended at $point with conditions $conds"))
    end

end

function GetNewPoint(point, line)

    # Finding the gradient and intercept of the perpendicular line from the point to the boundary
    perpGrad = -1/line.m
    perpIntercept = -perpGrad * point[1] + point[2]

    # Defining the perpendicular line
    perpLine = BoundLine(perpGrad, perpIntercept)

    # Getting the coordinates
    x = (line.c - perpLine.c) / (perpLine.m - line.m)
    y = LineEquation(x, line)

    point = [x,y]

end

function NearestLine(point, lines, a)

    # Finding the points at which the top and bottom triangle bounds intersect eachother

    # Top triangle vertex
    xi1 = (lines[2].c - lines[1].c) / (lines[1].m - lines[2].m)
    yi1 = LineEquation(xi1, lines[1])

    # Bottom triangle vertex
    xi2 = (lines[4].c - lines[3].c) / (lines[3].m - lines[4].m)
    yi2 = LineEquation(xi2, lines[3])
    
    # Returning the closest boundary line that the point is violating
    if point[2] > yi1
        # Checking that the new point will not be in breach of the constraints
        if (point[1] < xi1) & (GetNewPoint(point, lines[1])[1] > Constraints(a))
            return lines[1]
        else
            return lines[2]
        end
    elseif point[2] < yi2
        if (point[1] < xi2) & (LineEquation(point[1], lines[3]) > Constraints(a))
            return lines[3]
        else
            return lines[4]
        end
    end

end

function ReturnToConstraints(param, a)

    bound = Constraints(a)

    if param < bound
        param = bound
    end

    return param

end

