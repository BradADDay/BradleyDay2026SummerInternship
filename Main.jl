using Gradus
using Plots
using ColorSchemes
using SpectralFitting

include("LampPostModelFit.jl")
include("ParameterVariations.jl")
include("Defaults.jl")

# =======================================================================
# Parameter space
# =======================================================================

# Looping through the values for each variable and computing the line profile
#ParamLoop("θ", range(5, 85, 4), copy(defaultSetupDict); line=false, render=true, imageSize=(400,300))

# =======================================================================
# Fitting
# =======================================================================

# Model parameters
setupDict = Dict((
                  "M"   => 1., 
                  "a"   => 0.998, 
                  "α13" => 0., 
                  "α22" => 0., 
                  "α52" => 0.,
                  "ϵ3"  => 0., 
                  "θ"   => 70., 
                  "h"   => 20.))

# Generating bins
bins = collect(range(0.2, 1.5, 100))

# Simulating a noisy emission line
flux = JohannsenParamVar(setupDict, bins, ComputeLineProfile; 
                render = false, minrₑ = -1., maxrₑ = 400., 
                numrₑ = 100)

noise = 0.002

flux += rand(-noise:1e-8:noise, length(flux))
flux[flux.<0] .= 0

# Putting the flux into a data object
data = InjectiveData(bins, flux, name="Data")

# Plotting
plot(data, markersize=3)
display(plot!(xlabel="Energy", ylabel="Flux (arb. units)"))

# Setting up the fitting problem
model = LampPostJohannsen()
prob = FittingProblem(model => data)

# Fitting
println("+ Fitting...")
result = SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, maxIter = 25, verbose=true)

# Plotting
plot!(result)
