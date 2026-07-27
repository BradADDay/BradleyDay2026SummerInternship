using CSV, DataFrames, Interpolations, Gradus

const DeformationBoundsTable = DataFrame(CSV.File("Code/utils/DeformationBounds.csv"))

# ====================================================================================
# Lines
# ====================================================================================

"""
    Gradient(coord1, coord2)

Calculate the gradient of the line connecting two points
"""
function Gradient(coord1::AbstractArray, coord2::AbstractArray)::Real
    diff = coord1 .- coord2
    diff[2]/diff[1]
end

"""
    Intercept(grad, coord)

Calculate the Intercept of the line passing through a point with a given gradient
"""
Intercept(grad::Real, coord::AbstractArray)::Real = coord[2] - grad * coord[1]

"""
    BoundLine{T<:Real}

A structure to define a line in 2D
"""
struct BoundLine{T<:Real}
    "Gradient"
    m::T
    "Intercept"
    c::T
end

"""
Constructor from two coordinates
"""
function BoundLine(coord1::AbstractArray, coord2::AbstractArray)
    m = Gradient(coord1, coord2)
    c = Intercept(m, coord1)
    BoundLine(m, c)
end

"""
    LineEquation(x, line)

Get the y coordinate of a point on a line from a given x coordinate.
Uses the `BoundLine` struct.
"""
function LineEquation(x::Real, line::BoundLine)::Real
    line.m * x .+ line.c
end

"""
    Intersection(line1, line2)

Find the intersection point of two `BoundLine` objects
"""
function Intersection(line1::BoundLine, line2::BoundLine)::Real
    (line2.c - line1.c) / (line1.m - line2.m)
end

"""
    GetLineParams(a, df)

Get the Intercepts and Gradients for the boundary lines for a given a from the data table.
"""
function GetLineParams(a::Real, df::DataFrame=DeformationBoundsTable)

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

"""
    GetLines(a, df)

Get the line objects defining the parameter space boundaries for a given spin value.
"""
function GetLines(a::Real, df::DataFrame=DeformationBoundsTable)
    
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

"""
    ValidityCheckISCO(m)

Check if a given Johannsen metric is valid by checking for naked singularities, the 
presence of an ISCO, and the boundaries. If it is valid then the value of the ISCO is 
returned
"""
function ValidityCheckISCO(m::AbstractMetric)::Real
    
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

"""
    IsValid(ϵ3, α13, a)

Check if a parameter combination produces a valid Johannsen metric by explicitly checking
for naked singularities, an ISCO, or if it is outside of predefined bounds.

See also [`QuickIsValid`](@ref)
"""
function IsValid(ϵ3::Real, α13::Real, a::Real)::Bool
    m = JohannsenMetric(a=a, ϵ3=ϵ3, α13=α13)
    
    if is_naked_singularity(m) | is_no_isco(m) | (ϵ3 < Constraints(a)) | (α13 < Constraints(a)) | (abs(a) > 0.998)
        return false
    else 
        return true
    end

end

"""
    QuickIsValid(ϵ3, α13, a)

Check if a parameter combination produces a valid Johannsen metric by comparing against the 
defined boundary lines.

See also [`IsValid`](@ref)
"""
function QuickIsValid(ϵ3::Real, α13::Real, a::Real, df::DataFrame=DeformationBoundsTable; 
        sgns::AbstractArray=[1,1,-1,-1,-1,-1]
    )::Bool

    if (ϵ3 < Constraints(a)) | (α13 < Constraints(a)) | (abs(a) > 0.998)
        return false
    end
    
    conds = GetValidityConditions(ϵ3, α13, a, df; sgns=sgns)
    
    # Checking if the parameter combination is within any of the 4 forbidden regions
    if conds[1] | conds[2] | conds[3] | conds[4]
        return false
    else
        return true
    end

end

"""
    GetValidityConditions(ϵ3, α13, a)

Check a parameter combination against the pre defined boundary lines to see if it is valid.
"""
function GetValidityConditions(
        ϵ3::Real, α13::Real, a::Real, 
        df::DataFrame=DeformationBoundsTable; sgns::AbstractArray=[1,1,-1,-1,-1,-1]
    )

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

"""
    is_no_isco(m)

Check if there is a valid ISCO for the input metric
"""
function is_no_isco(m::AbstractMetric)::Bool
    try 
        Gradus.isco(m)
        return false
    catch
        return true
    end
end

"""
    is_no_isco(ϵ3, α13, a)

Check if there is a valid ISCO for the input parameter combination
"""
function is_no_isco(ϵ3::Real, α13::Real, a::Real)::Bool
    is_no_isco(JohannsenMetric(ϵ3=ϵ3, α13=α13, a=a))
end

@doc raw"""
    Constraints(a)

The boundary line for the Johannsen deformation parameters α13 and ϵ3.
```math
-(1+\sqrt{1 - a^2)})^3
```
"""
Constraints(a::Real)::Real = -(1+sqrt(1 - a^2))^3

"""
    IsValidFit(fit)

Check is a `SpectralFitting` Johannsen metric `FitResult` is valid 
"""
function IsValidFit(fit)::Bool
    # Check if a SpectralFitting FitResult is valid
    u=fit.u
    IsValid(u[7], u[6], u[3])
end

# ====================================================================================
# Corrections
# ====================================================================================

"""
    FindNearestSafePoint(point, a; maxIter=1000)

Find the nearest safe point in the Johannsen parameter space by finding the nearest point 
on the offending boundary lines
"""
function FindNearestSafePoint(
        point::AbstractArray, a::Real; maxIter::Int=1000
    )::AbstractArray
    
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

"""
    GetNewPoint(point, line)

Find the closest point from a point to a line in 2D.
"""
function GetNewPoint(point::AbstractArray, line::BoundLine)::AbstractArray

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

"""
    NearestLine(point, lines, a)

Find the nearest line from a given point out a selection of lines. The output will always 
lead [`GetNewPoint`](@ref) to return a valid point in the parameter space.
"""
function NearestLine(
        point::AbstractArray, lines::AbstractArray{BoundLine}, a::Real
    )::BoundLine

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

"""
    ReturnToConstraints(param, a)

If a parameter is out of bounds for a given spin value, returning it to the nearest valid 
value (on the boundary line).
"""
function ReturnToConstraints(param::Real, a::Real)::Real

    bound = Constraints(a)

    if param < bound
        param = bound
    end

    param

end

