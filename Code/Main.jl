## ======================================================================

using Gradus
using Plots
using SpectralFitting
using Colors
using Dates

pyplot()

include("FittingModels.jl")
include("ParameterVariations.jl")
include("Defaults.jl")

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    gridalpha=0.,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9))
)

E = 6.4
a = 0.998
h = 11.
θ = 45.
α13 = 5.
ϵ3 = -1.

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

bins = range(3, 10, 1000)

model = XS_LampPostJohannsen(;K=1., E=E, a=a, h=h, θ=θ, α13=α13, ϵ3=ϵ3)
flux = invokemodel!(bins, model)
plot(bins[1:end-1], flux/maximum(flux); label="Table", xlabel="Energy (keV)", ylabel="(arb units)")

flux1 = JohannsenParamVar(setupDict, bins/6.4, ComputeLineProfile)
plot!(bins, flux1/maximum(flux1); label="Gradus")

##
for α13 in 0:0.2:1., ϵ3 in 0:0.2:1.
     α13 = -α13
     ϵ3 = -ϵ3
     println(α13, " ", ϵ3)
     try
          # Model parameters
          setupDict = Dict((
                         "M"   => 1., 
                         "a"   => a, 
                         "α13" => α13, 
                         "α22" => 0., 
                         "α52" => 0.,
                         "ϵ3"  => ϵ3, 
                         "θ"   => θ, 
                         "h"   => h))


          bins = range(1, 9, 999)

          flux1 = JohannsenParamVar(setupDict, bins/6.4, ComputeLineProfile)

     catch err
          println("Failed")
          open("./err.txt", "a") do io
               write(io, "$α13, $ϵ3\n$err\n")
          end
     end
end
