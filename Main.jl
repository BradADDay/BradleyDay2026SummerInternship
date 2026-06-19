## ======================================================================

using Gradus
using Plots
using SpectralFitting
using Colors
using Dates

gr()

include("LampPostModelFit.jl")
include("ParameterVariations.jl")
include("Defaults.jl")

# =======================================================================
# Parameter space
# =======================================================================

# Looping through the values for each variable and computing the line profile
# ParamLoop("a", range(0, 0.95, 4), copy(defaultSetupDict); line=false, render=true, imageSize=(40,30))

## =======================================================================
# Fitting
# =======================================================================

t0 = now()

# Model parameters
setupDict = Dict((
                  "M"   => 1., 
                  "a"   => 0.8, 
                  "α13" => 20., 
                  "α22" => 0., 
                  "α52" => 0.,
                  "ϵ3"  => 4., 
                  "θ"   => 70., 
                  "h"   => 20.))

# Generating bins
bins = collect(range(0.2, 1.5, 100))

# Simulating a noisy emission line
flux = JohannsenParamVar(setupDict, bins, ComputeLineProfile; 
                render = false, minrₑ = -1., maxrₑ = 400., 
                numrₑ = 100)

scatter(bins, flux, label="True", minorticks=4, gridalpha=0.5, 
     minorgrid=true, minorgridalpha=0.3, xlabel="Energy", 
     ylabel="Flux (arb. units)", markersize=2)

noise = 0.002
noisyFlux = flux + rand(-noise:1e-8:noise, length(flux))
noisyFlux[noisyFlux.<0] .= 0

# Putting the flux into a data object
data = InjectiveData(bins, noisyFlux, name="Noisy")

# Plotting
display(plot!(data, markersize=3))

# Setting up the fitting problem
model = LampPostJohannsen()
prob = FittingProblem(model => data)

# Fitting
println("+ Fitting...")
result = SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true)

# Plotting
plot!(result)

println("Time elapsed: $(now() - t0)")