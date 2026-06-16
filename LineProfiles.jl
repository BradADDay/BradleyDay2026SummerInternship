using Gradus
using Plots
using ColorSchemes

gr()

println(Threads.nthreads(:default))

# =======================================================================
# Functions
# =======================================================================

function locate(value, array)
    findall(x -> x==value, array)[1]
end

function generateConfig(variables, setupDict)
    """
    Generate a string denoting the configuration of a given data run.
    """
    config = ""

    # Looping through the variables
    for i in variables
        config *= "$i = $(round(setupDict[i]; digits= 2))"
        if i != variables[end]
            config *= ", "
        end
    end
    return config
end

function computeLineProfile(m, x, height, bins = range(0.0, 1.5, 180))

    # Disk
    d = ThinDisc(0.0, Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    bins, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
            method=TransferFunctionMethod())

    return bins, flux
end

function renderImage(m, x, λ_max, imageSize=[40,30])

    # Disk
    d = ThinDisc(0.0, 15.0)

    # Redshift point function
    redshift = ConstPointFunctions.redshift(m, x)
    redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()
    
    # Rendering the image
    α, β, image = rendergeodesics(
            m, x, d, λ_max, pf = redshiftGeometry,
            # image parameters
            image_width = imageSize[1], image_height = imageSize[2],
            αlims = (-20, 20), βlims = (-15, 15), verbose = true)

    return α, β, image
end

function paramVar(setup, bins, config, line, render)

    returns = []

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(setup["incl"]), 0.0)
    λ_max = 2x[2]  # Max affine time ~2x[2]

    # Instantiating the metric
    m = JohannsenMetric(setup["M"], setup["a"], setup["α13"], setup["α22"], setup["α52"], setup["ϵ3"])

    # Rendering the image
    if render
        println("+ Rendering for $config")
        α, β, image = renderImage(m, x, λ_max, [200, 150])

        append!(returns, [α, β, image])
    end
    
    # Computing the line profile
    if line
        println("+ Computing line profile for $config")
        bins, flux = computeLineProfile(m, x, setup["h"], bins)
        append!(returns, [bins, flux])
    end

    return returns

end

function paramLoop(parameter, values, setupDict; line=true, render=true)

    # Setup arrays for storage
    configs = []
    binVars = []
    fluxes = []
    αs = []
    βs = []
    images = []

    # Initialising bins
    bins = collect(range(0.1, 1.5, 180))

    for i in values
        try
            setupDict[parameter] = i

            config = generateConfig([parameter], setupDict)
            result = paramVar(setupDict, bins, config, line, render)
            
            if line & render
                bins, flux, α, β, image = result
                append!(binVars, [bins])
                append!(fluxes, [flux])
                append!(αs, [α])
                append!(βs, [β])
                append!(images, [image])
            elseif line
                bins, flux = result
                append!(binVars, [bins])
                append!(fluxes, [flux])
            elseif render
                α, β, image = result
                append!(αs, [α])
                append!(βs, [β])
                append!(images, [image])
            end

            append!(configs, [config])
        catch
            println("Value of $i failed to compute")
        end
    end

    if line & render
        return configs, binVars, fluxes, αs, βs, images
    elseif line
        return configs, binVars, fluxes
    elseif render
        return configs, αs, βs, images
    end
end

# =======================================================================
# Setup
# =======================================================================

function defaultParameters()
    # Metric Parameters
    M = 1.0     # Mass
    a = 0.998   # Spin
    # Perturbations
    α13 = 0.0
    α22 = 0.0
    α52 = 0.0
    ϵ3  = 0.0

    # BH parameters
    incl = 60.0 # Inclination, degrees
    h    = 10.0 # Corona height

    # Dictionary for easy access and modification of the model parameters
    setupDict = Dict("incl" => incl, 
                    "α13"  => α13, 
                    "M"    => M, 
                    "α22"  => α22, 
                    "ϵ3"   => ϵ3, 
                    "a"    => a, 
                    "h"    => h, 
                    "α52"  => α52)

    return setupDict
end

# =======================================================================
# Function calls
# =======================================================================

# Variables to be investigated and the ranges to do so within
variables = ["incl", "α13", "ϵ3", "a", "h"]
# ranges = [5:10:85, 0:10:50, 0:10:50, 0.05:0.1:0.95, 5:5:30]
ranges = [range(5, 85, 4), range(0, 50, 4), range(0, 30, 4), range(0.05, 0.95, 4), range(5, 30, 4)]

#= variables = ["α22"]
ranges = [0:10:50] =#

# Looping through the variables
for i in 1:length(ranges)

    setupDict = defaultParameters()

    range = ranges[i]
    variable = variables[i]

    # Looping through the values for each variable and computing the line profile
    # configs, binVars, fluxes = paramLoop(variable, range, setupDict; line=false)
    configs, αs, βs, images = paramLoop(variable, range, setupDict; line=false)

    display(heatmap(αs, βs, images, aspect_ratio = 1; layout = length(αs), title=[i for j in 1:1, i in configs], titleloc=:center))

    #= plot(xlabel="g", ylabel="Flux (Arbitrary)", title="Variations of $variable")

    for i in 1:length(binVars)
        display(plot!(binVars[i], fluxes[i]; label=configs[i], palette=:tab10, lw=3))
        println("($i) $(configs[i])")
    end =#
end
