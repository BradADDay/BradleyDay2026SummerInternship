using Gradus
using Plots
using ColorSchemes
using SpectralFitting
using XSPECModels


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

function computeLineProfile(m, x, height, bins = range(0.0, 1.5, 180))

    # Disk
    d = ThinDisc(10., Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
            method=TransferFunctionMethod(), minrₑ=10.)

    #= profile(r) = r^-2

    _, flux = lineprofile(bins, profile, m, x, d; verbose=true,
            method=TransferFunctionMethod(), minrₑ=10.) =#

    return flux, Gradus.isco(m), profile
end

function renderImage(m, x, λ_max, imageSize=(40,30))

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

    return α, β, image
end

function paramVar(setup, bins, config; line=true, render=true)

    returns = Dict()

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(setup["incl"]), 0.0)
    λ_max = 2x[2]  # Max affine time ~2x[2]

    # Instantiating the metric
    m = JohannsenMetric(setup["M"], setup["a"], setup["α13"], setup["α22"], setup["α52"], setup["ϵ3"])

    # Rendering the image
    if render
        println("+ Rendering for $config")
        α, β, image = renderImage(m, x, λ_max, (400, 300))
        returns["α"] = α
        returns["β"] = β
        returns["image"] = image
    end
    
    # Computing the line profile
    if line
        println("+ Computing line profile for $config")
        flux, isco, profile = computeLineProfile(m, x, setup["h"], bins)
        returns["flux"] = flux
        returns["profile"] = profile
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
            # Base.show_backtrace(stdout, backtrace())
        end
    end

    plotting(parameter, results, bins; line, render)
end

function plotting(parameter, results, bins; line, render)
    configs = results["config"]

    if render
        # Plotting the heatmaps as subplots
        display(heatmap(results["α"], results["β"], results["image"], aspect_ratio = 1; 
                layout = length(configs), title=[i for j in 1:1, i in configs], titleloc=:center))

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
const a = 0.   # Spin

# Perturbation Defaults
const α13 = 0.
const α22 = 0.
const α52 = 0.
const ϵ3  = 0.

# BH Defaults
const incl = 60. # Inclination, degrees
const h    = 100. # Corona height

# Dictionary for easy access and modification of the model parameters
defaultSetupDict = Dict((["incl", incl], ["α13", α13], ["M", M], ["α22", α22], ["ϵ3", ϵ3], ["a", a], ["h", h], ["α52", α52]))

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

bins = collect(range(0.5, 1.5, 180))
result = paramVar(defaultSetupDict, bins, "Default Parameters"; render=false)
flux = result["flux"]
profile = result["profile"]

data = InjectiveData(bins, flux, name="Data")

model = XS_DiskLine(θ=FitParam(60.0), lineE=FitParam(1.), inner_r=FitParam(10., frozen=true), outer_r=FitParam(50.0, frozen=false))

prob = FittingProblem(model => data)
result = SpectralFitting.fit(prob, LevenbergMarquadt())

plot(data, markersize=2)
plot!(result)