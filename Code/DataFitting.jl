using BSON: @save, @load
using BenchmarkTools

# Loading in files
include("utils/FittingUtils.jl")

# @load "output/testResultXRISM.bson" result

# @model function mcmc_model(objective, stddev, f)
#     K ~ Normal(20.0, 1.0)
#     a ~ Normal(2.2, 0.3)
#     ηH ~ truncated(Normal(0.5, 0.1); lower = 0)
#     pred = f(K, a, ηH)
#     return objective ~ MvNormal(pred, stddev)
# end

# ==========================================================================================
# XRISM
# ==========================================================================================

# Loading the data files
data = LoadData(
    "data/xa000125000xtd_src_rebinned.pi", 
    "data/xa000125000xtd_bgd_rebinned.pi", 
    "data/xa000125000xtd_p031100010_src.rmf", 
    "data/xa000125000xtd_p031100010_ptsrc.arf";
    dataRange=(3,10)
)

# TODO: Maybe try changing the bin widths if this doesnt work?

# Lamp post line profile model
# TODO: Add R_OUT as a parameter and increase range of possible corona heights?
Γ = 3.5
θ = 15.

LP = XS_LampPostJohannsen(;
    E = FitParam(1., lower_limit=1., upper_limit=10, frozen=true),
    a = FitParam(0.998, lower_limit=0, upper_limit=0.998, frozen=false),
    h = FitParam(19., lower_limit=3., upper_limit=19., frozen=false),
    θ = FitParam(θ, lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(0.01, lower_limit = 0., upper_limit=Inf, frozen=false),
    Γ = FitParam(Γ, lower_limit = 1., upper_limit = 5., frozen = false),
    A_Fe = FitParam(1., lower_limit = 0., upper_limit = 2., frozen = true),
    logXi = FitParam(0.5, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(30., lower_limit=1., upper_limit=30., frozen = true), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(10., lower_limit=0.0, upper_limit=30.0, frozen=false)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(4e4, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=0., upper_limit=5., frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LP) 
#PL + 
# Defining the model
model = AB * (PL + ConvModel(SP))

# Defining the problem and binding parameters
prob = FittingProblem(model => data)
bind!(prob, (1, :a1, :a) => (1, :a2, :Γ))
bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination))

# Fitting
result = SpectralFitting.fit(
    prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1e3)
)

display(result)
    
plt = PlotResult(data, result[1])
display(plt)

close("all")

@save "output/FirstFullResult.bson" result

## =========================================================================================
# NUStar
# ==========================================================================================

# List of available datasets
files = [
    "nu80402315002", 
    "nu80402315004", 
    "nu80402315006", 
    "nu80402315008", 
    "nu80402315010", 
    "nu80402315012", 
    "nu80502304002", 
    "nu80502304004", 
    "nu80502304006"
]

_, johannsenResult = FitNUStar(files[1]; kerr=false)

@save "output/testfit2.bson" johannsenResult

close("all")

##
# Testing

Γ = 3.5
θ = 15.

LP = XS_LampPostJohannsen(;
    E = FitParam(1., lower_limit=1., upper_limit=10, frozen=true),
    a = FitParam(0.998, lower_limit=0, upper_limit=0.998, frozen=false),
    h = FitParam(19., lower_limit=3., upper_limit=19., frozen=false),
    θ = FitParam(θ, lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(0.01, lower_limit = 0., upper_limit=Inf, frozen=false),
    Γ = FitParam(Γ, lower_limit = 1., upper_limit = 5., frozen = false),
    A_Fe = FitParam(1., lower_limit = 0., upper_limit = 2., frozen = true),
    logXi = FitParam(0.5, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(30., lower_limit=1., upper_limit=30., frozen = true), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(10., lower_limit=0.0, upper_limit=3.0, frozen=true)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(4e4, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=0., upper_limit=2., frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LP) 
#PL + 
# Defining the model
model = AB * (PL + ConvModel(SP))

bins = 3:0.01:10

ConvSPFlux = invokemodel(bins, ConvModel(SP))
SPFlux = invokemodel(bins, SP)
ABFlux = invokemodel(bins, AB)
PLFlux = invokemodel(bins, PL)

flux = ABFlux .* (SPFlux .+ PLFlux)
convflux = ABFlux .* (ConvSPFlux .+ PLFlux)

plot(bins, flux)
plot!(bins, convflux)


