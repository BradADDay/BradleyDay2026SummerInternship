## ======================================================================
using Dates
using SpectralFitting

include("utils/FittingModels.jl")
include("utils/PlottingDefaults.jl")
include("utils/ParameterVariations.jl")

E = 6.4
a = 0.9
h = 12.
θ = 76.
α13 = 0.
ϵ3 = 0.
K = 30.

bins = range(3, 10, 1000)

# Position of the observer
x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

# Instantiating the metric
m = JohannsenMetric(M=1., a=a, α13=α13, ϵ3=ϵ3)
# print("ISCO: $(Gradus.isco(m))")

# Setting the inner radius to the ISCO if the entered value is <0
minrₑ = Gradus.isco(m)

# Disk
d = ThinDisc(0., Inf)

# Setting up the model and emissivity profile
model = LampPostModel(h = h)
profile = emissivity_profile(m, d, model)

# Computing the line profile
_, flux = lineprofile(m, x, d, profile; verbose=false, bins=bins/E, 
          method=TransferFunctionMethod())

plot(bins, flux/maximum(flux))

model = XS_LampPostJohannsen(;K=K, E=E, a=a, h=h, θ=θ, α13=α13, ϵ3=ϵ3)
flux2 = invokemodel!(bins, model)

println(IsValid(ϵ3, α13, a))
println(flux2[1:20])

plot!(bins[1:end-1], flux2/maximum(flux2))

##

E = 6.4
a = 0.998
h = 10.
θ = 15.
α13 = 0.
ϵ3 = 0.
K = 10.

model = XS_LampPostJohannsen(;K=K, E=E, a=a, h=h, θ=θ, α13=α13, ϵ3=ϵ3)
flux = invokemodel!(bins, model)
flux2 = invokemodel(bins, PowerLaw(K=FitParam(100.), a=FitParam(0.5)))

plot(bins[1:end-1], flux; label="Table", xlabel="Energy (keV)", ylabel="(arb units)")
plot!(bins[1:end-1], flux2)
plot!(bins[1:end-1], flux2 .+ flux)
display(vline!([6.4]))

GC.gc()
close("all")


##

using Dates
using SpectralFitting

include("utils/FittingModels.jl")
include("utils/PlottingDefaults.jl")
include("utils/ParameterVariations.jl")


function GetLineProfile(bins, a, h, θ, α13, ϵ3)

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

    # Instantiating the metric
    m = JohannsenMetric(M=1., a=a, α13=α13, ϵ3=ϵ3)

    # Disk
    d = ThinDisc(0., Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = h)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=false, bins=bins, 
            method=TransferFunctionMethod()
    )

    return flux

end

bins = range(3, 10, 1000)

E = 6.4
a, h, θ, α13, ϵ3 = (0.998, 3.0, 5.0, -0.7018823703857276, -1.2018823703857278)
K = 10.


setupDict = Dict((
    "M"   => 1., 
    "a"   => a, 
    "α13" => α13, 
    "α22" => 0., 
    "α52" => 0.,
    "ϵ3"  => ϵ3, 
    "θ"   => θ, 
    "h"   => h
))

flux2 = GetLineProfile(bins, a, h, θ, α13, ϵ3)

plot(bins[1:end-1], flux2/maximum(flux2))

ϵ3 +=1
m = JohannsenMetric(M=1., a=a, ϵ3=ϵ3, α13=α13)
is_no_isco(m)
is_naked_singularity(m)