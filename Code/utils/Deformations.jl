using CSV, DataFrames, Interpolations, Gradus

# export GetLines, ValidityCheckISCO, IsValid, QuickIsValid, is_no_isco, Constraints, Gradient, Intercept, BoundLine, LineEquation

const DeformationBoundsTable = DataFrame(CSV.File("Code/utils/DeformationBounds.csv"))

function GetLineParams(a, df=DeformationBoundsTable)

    # Read the file and separate the spins
    df = copy(df)
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

function ValidityCheckISCO(m)
    """
    Checking if the metric is valid
        returns the ISCO if so
        returns an error value if not
    """
    if is_naked_singularity(m)
        # Naked singularity
        return -3
    elseif ((m.α13 < Constraints(m.a)) | (m.ϵ3 < Constraints(m.a)))
        # Abnormal exterior region
        return -2
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
    
    if is_naked_singularity(m) | is_no_isco(m) | (ϵ3 < Constraints(a)) | (α13 < Constraints(a))# | (abs(a) > 0.998)
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

    # Getting the boundary lines and defining conditions
    lines = GetLines(a, df)
    conds = fill(true, length(lines))

    # Checking if the parameter combination is within the forbidden region for every line
    for i in eachindex(lines)
        if !isnothing(lines[i])
            conds[i] = sgns[i]*α13 > sgns[i]*(LineEquation(ϵ3, lines[i]))
        end
    end
    
    # Checking if the parameter combination is within any of the 4 forbidden regions
    if (conds[1] & conds[2]) | (conds[3] & conds[4]) | conds[5] | conds[6]
        return false
    else
        return true
    end

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

function is_no_isco(ϵ3, α13, a)
    
    is_no_isco(JohannsenMetric(ϵ3=ϵ3, α13=α13, a=a))
    
end

# The lower limit for the deformation parameters α13, ϵ3 for a given spin
Constraints(a) = -(1+sqrt(1 - a^2))^3

# Defining line parameters
function Gradient(coord1, coord2)
    diff = coord1 .- coord2
    diff[2]/diff[1]
end

function Intercept(grad, coord)
    coord[2] - grad * coord[1]
end

# Line struct
struct BoundLine{T<:Real}
    "Gradient"
    m::T
    "Intercept"
    c::T
end

# Constructor
function BoundLine(coord1::AbstractArray, coord2::AbstractArray)
    m = Gradient(coord1, coord2)
    c = Intercept(m, coord1)
    BoundLine(m, c)
end

# Get a point on the line
function LineEquation(x, line::BoundLine)
    line.m * x .+ line.c
end

function IsValidFit(fit)
    u=fit.u
    IsValid(u[7], u[6], u[3])
end