using Gradus
using Plots
using ColorSchemes
using SpectralFitting

include("LampPostModelFit.jl")

# =======================================================================
# Functions
# =======================================================================

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

function computeLineProfile(m, x, height, bins = range(0.0, 1.5, 180); minrₑ=-1., maxrₑ=400., numrₑ=100, kwargs...)

    if minrₑ < 0.
        minrₑ = Gradus.isco(m)
    end

    # Disk
    d = ThinDisc(minrₑ, Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
            method=TransferFunctionMethod(), minrₑ=minrₑ, maxrₑ=maxrₑ, numrₑ=200)

    return flux
end

function renderImage(m, x, λ_max; imageSize=(40,30), kwargs...)

    # Disk
    d = ThinDisc(0.0, 15.0)

    # Redshift point function
    redshift = ConstPointFunctions.redshift(m, x)
    redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()

    time_coord = PointFunction((m, gp, λ) -> gp.x[1])
    pfGeometry = time_coord ∘ ConstPointFunctions.filter_intersected()

    # Rendering the image
    α, β, image = rendergeodesics(
            m, x, d, λ_max, pf = redshiftGeometry,
            # image parameters
            image_width = imageSize[1], image_height = imageSize[2],
            αlims = (-20, 20), βlims = (-15, 15), verbose = true)
    
    display(heatmap(α, β, image, aspect_ratio = :equal; clims=(0., 1.5), cmap=:redsblues, ))

    return α, β, image
end

function paramVar(setup, bins, config; render=true, line=true, kwargs...)

    returns = Dict()

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(setup["θ"]), 0.0)
    λ_max = 2x[2]  # Max affine time ~2x[2]

    # Instantiating the metric
    m = JohannsenMetric(setup["M"], setup["a"], setup["α13"], setup["α22"], setup["α52"], setup["ϵ3"])

    # Rendering the image
    if render
        println("+ Rendering for $config")
        α, β, image = renderImage(m, x, λ_max; kwargs...)
        returns["α"] = α
        returns["β"] = β
        returns["image"] = image
    end
    
    # Computing the line profile
    if line
        println("+ Computing line profile for $config")
        flux = computeLineProfile(m, x, setup["h"], bins; kwargs...)
        returns["flux"] = flux
    end

    return returns

end

function paramLoop(parameter, values, setupDict, bins; 
                   line=true, render=true)

    # Setup arrays for storage
    results = Dict("flux" => [],
                   "α" => [],
                   "β" => [],
                   "image" => [],
                   "config" => [])

    keys = []
    if line
        append!(keys, ["flux"])
    end
    if render
        append!(keys, ["α", "β", "image"])
    end

    for i in values
        try
            setupDict[parameter] = i

            config = generateConfig([parameter], setupDict)
            result = paramVar(setupDict, bins, config; line, render)

            for key in keys
                append!(results[key], [result[key]])
            end

            append!(results["config"], [config])
        catch err
            println("Value of $i failed to compute")
            println(err)
        end
    end

    plotting(parameter, results, bins; line, render)
end

function plotting(parameter, results, bins; line, render)
    configs = results["config"]

    if render
        # Plotting the heatmaps as subplots
        display(heatmap(results["α"], results["β"], results["image"], aspect_ratio = 1; 
                layout = length(configs), title=[i for j in 1:1, i in configs], titleloc=:center, cmap=:redsblues))

        # Plotting the heatmaps individually
        for i in eachindex(configs)
            display(heatmap(results["α"][i], results["β"][i], results["image"][i], aspect_ratio = :equal, title=configs[i]; clims=(0., 1.5)))
            savefig("testimages/$(configs[i]).png")
        end
    end

    if line
        # Plotting the line profile
        plot(xlabel="g", ylabel="Flux (Arbitrary)")

        for i in eachindex(configs)
            plot!(bins, results["flux"][i]; label=configs[i], palette=:tab10, lw=3)
        end

        display(plot!(title="Variations of $parameter"))
    end
end

# =======================================================================
# Setup
# =======================================================================

# Metric Defaults
const M = 1.0     # Mass
const a = 0.998   # Spin

# Perturbation Defaults
const α13 = 0.
const α22 = 0.
const α52 = 0.
const ϵ3  = 0.

# BH Defaults
const θ = 60. # Inclination, degrees
const h = 9. # Corona height

# Dictionary for easy access and modification of the model parameters
defaultSetupDict = Dict((["θ", θ], ["α13", α13], ["M", M], ["α22", α22], ["ϵ3", ϵ3], ["a", a], ["h", h], ["α52", α52]))

# Parameter Variations
variables = ["ϵ3"]
values = [32.2:0.1:34.8]

# =======================================================================
# Function calls
# =======================================================================

bins = collect(range(0.1, 1.5, 180))

#= # Looping through the variables
for i in eachindex(values)
    println()

    setupDict = copy(defaultSetupDict)

    value = values[i]
    variable = variables[i]

    # Looping through the values for each variable and computing the line profile
    paramLoop(variable, value, setupDict, bins; line=true, render=true)
end =#

bins = collect(range(0.2, 1.5, 100))
flux = paramVar(defaultSetupDict, bins, "Default Parameters"; 
                render = false, minrₑ = -1., maxrₑ = 400., 
                numrₑ = 10)["flux"]

flux += rand(-0.001:1e-7:0.001, length(flux))
flux[flux.<0] .= 0

data = InjectiveData(bins, flux, name="Data")

model = LampPostJohannsen()

prob = FittingProblem(model => data)

println("+ Fitting...")
result = SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, maxIter = 15, verbose=true)

plot(data, markersize=3)
plot!(result)
plot!(xlabel="Energy", ylabel="Flux (arb. units)")

flux = paramVar(defaultSetupDict, bins, "Default Parameters"; 
                render = false, minrₑ = -1., maxrₑ = 400., 
                numrₑ = 10)["flux"]