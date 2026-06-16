using Gradus
using Plots

Plots.default(show = true)

gr()

# =======================================================================
# Functions
# =======================================================================

function computeLineProfile(m, x, d, height, bins = range(0.0, 1.5, 180))

    # Setting up the model and profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    bins, flux = lineprofile(m, x, d, profile, bins=bins, 
            method=TransferFunctionMethod())

    return bins, flux
end

function renderImage(m, x, d, λ_max, imageSize=[40,30])

    # Disk
    d = ThinDisc(0.0, 15)

    # Redshift point function
    redshift = ConstPointFunctions.redshift(m, x)
    redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()
    
    α, β, image = rendergeodesics(
            m, x, d, λ_max, pf = redshiftGeometry,
            # image parameters
            image_width = imageSize[1], image_height = imageSize[2],
            αlims = (-20, 20), βlims = (-15, 15), verbose = true
            )

    display(heatmap(α, β, image, aspect_ratio = 1))
end

# =======================================================================
# Setup
# =======================================================================

# Metric Parameters
M = 1.0     # Mass
a = 0.998M  # Spin (Fraction of M)
# Perturbations
α13 = 0.0
α22 = 0.0
α52 = 0.0
ϵ3 = 0.0

# BH parameters
i = deg2rad(70) # Inclination
h = 10.0        # Corona height

# Disk
d = ThinDisc(0.0, Inf)

# =======================================================================
# Function calls
# =======================================================================

# Initialising the plot

binVars = []
fluxVars = []
labels = []

bin = collect(range(0.1, 1.5, 90))

for ϵ3 in collect(0.0:5.0:50.0)

    try
        # Position of the observer
        # (x[1], Distance, Inclination, x[4])
        x = SVector(0.0, 10000.0, i, 0.0)
        λ_max = 2x[2]  # Max affine time ~2x[2]

        # Instantiating the metric
        m = JohannsenMetric(M, a, α13, α22, α52, ϵ3)
        
        #= println("Computing line profile for ϵ3 = $ϵ3")
        bins, flux = computeLineProfile(m, x, d, h, bin)

        append!(binVars, [bins])
        append!(fluxVars, [flux])
        append!(labels, ϵ3) =#

        println("Rendering for ϵ3=$ϵ3")
        renderImage(m, x, d, λ_max, [200, 150])
    catch
        println("Failed!")
    end
end

#= plot(xlabel="g", ylabel="Flux (Arbitrary)")

for i in 1:length(binVars)
    display(plot!(binVars[i], fluxVars[i]; label="ϵ3=$(labels[i])"))
end =#
