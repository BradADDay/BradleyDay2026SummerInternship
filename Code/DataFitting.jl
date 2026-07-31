using BSON: @save, @load
using BenchmarkTools
using Dates

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

# Lamp post line profile model
# TODO: Increase range of possible corona heights?
Γ = 3.4767
θ = 6.

LP = XS_LampPostJohannsen(;
    E = FitParam(1.0255, lower_limit=1., upper_limit=1.05, frozen=false),
    a = FitParam(0.998, lower_limit=0, upper_limit=0.998, frozen=false),
    h = FitParam(61.145, lower_limit=3., upper_limit=100., frozen=false),
    θ = FitParam(θ, lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(1.1493e-7, lower_limit = 0., upper_limit=1, frozen=false),
    Γ = FitParam(Γ, lower_limit = 1., upper_limit = 5., frozen = false),
    A_Fe = FitParam(0.96318, lower_limit = 0., upper_limit = 20., frozen = false),
    logXi = FitParam(0.46734, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(19.928, lower_limit=1., upper_limit=30., frozen = false), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(9.5425, lower_limit=5.0, upper_limit=13.0, frozen=false)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(0.27368, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=0., upper_limit=5., frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LP) 

# Defining the model
model = AB * (PL + ConvModel(SP))

# Defining the problem and binding parameters
prob = FittingProblem(model => data)
bind!(prob, (1, :a1, :a) => (1, :a2, :Γ))
bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination))

# Fitting
result = SpectralFitting.fit(
    prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1000)
)

display(result)

plt = PlotResult(data, result[1])
display(plt)
savefig(plt, "output/FitResult_$(data.user_data.object)_$(now()).png")

close("all")

@save "output/FitResult_$(data.user_data.object)_$(now()).bson" result

##
# Testing

Γ = 3.463
θ = 5.

LP = XS_LampPostJohannsen(;
    E = FitParam(1.02, lower_limit=1., upper_limit=1.02, frozen=true),
    a = FitParam(0.998, lower_limit=0, upper_limit=0.998, frozen=false),
    h = FitParam(40., lower_limit=3., upper_limit=40., frozen=false),
    θ = FitParam(θ, lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(1.3664e-7, lower_limit = 0., upper_limit=1, frozen=false),
    Γ = FitParam(Γ, lower_limit = 1., upper_limit = 5., frozen = false),
    A_Fe = FitParam(5., lower_limit = 0., upper_limit = 20., frozen = false),
    logXi = FitParam(.5, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(20., lower_limit=1., upper_limit=30., frozen = false), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(10.545, lower_limit=5.0, upper_limit=13.0, frozen=false)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(0.20956, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=0., upper_limit=2., frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LP) 
#PL + 
# Defining the model
model = AB * (PL + ConvModel(SP))

bins = 3:0.01:10

SPFlux = invokemodel(bins, ConvModel(SP))
ABFlux = invokemodel(bins, AB)
PLFlux = invokemodel(bins, PL)

flux = 50000 * ABFlux .* (SPFlux .+ PLFlux)

plot(data)
plot!(bins, flux)


