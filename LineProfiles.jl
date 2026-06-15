using Gradus
using Plots

Plots.default(show = true)

gr()

function computeLineProfile(m, i, height, bins = range(0.0, 1.5, 180))
    # Position of the observer
    # (x[1], Distance, Inclination, x[4])
    x = SVector(0.0, 10000.0, i, 0.0)

    # Disk
    d = ThinDisc(0.0, 100.0)

    # Setting up the model and profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    bins, flux = lineprofile(m, x, d, profile; bins=bins, verbose=true)

    return bins, flux
end

# Metric Parameters
M = 1.0
a = 0.998
α13 = 1.0
α22 = 0.0
α52 = 0.0
ϵ3 = 0.0
i = deg2rad(45)
h = 10.0

# Instantiating the metric
m = JohannsenMetric(M, a, α13, α22, α52, ϵ3)

# Initialising the plot
plot(xlabel="g", ylabel="Flux (Arbitrary)")

hs = [1.5, 3.0, 5.0, 10.0, 15.0, 20.0]

for h in 10.0:1.0:11.0
    bins, flux = computeLineProfile(m, i, h)
    
    display(plot!(bins, flux; label="h = $h"))
end

