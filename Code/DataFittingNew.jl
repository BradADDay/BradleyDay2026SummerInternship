using BSON: @save, @load
using BenchmarkTools
using Dates
using JSON3

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

# Line profile
E = 6.4
a = 0.998
h = 70.
θ = 5.
EConv = 1.0272

# Spectrum
Kₛ = 9.9329e-8
Γ = 1.5
A_Fe = 0.
logXi = 0.45761
density = 19.888

# Photoelectric absorption
ηH = 9.6071

# Power Law
Kₚ = 0.27368

vars = Dict(
    "E"=>E,
    "a"=>a,
    "h"=>h,
    "θ"=>θ,
    "EConv"=>EConv,
    "Kₛ"=>Kₛ,
    "Γ"=>Γ,
    "A_Fe"=>A_Fe,
    "logXi"=>logXi,
    "density"=>density,
    "ηH"=>ηH,
    "Kₚ"=>Kₚ
)

LP = XS_LampPostJohannsen(;
    K = FitParam(1., lower_limit=0., upper_limit=Inf, frozen=false),
    E = FitParam(E, lower_limit=6.4, upper_limit=7., frozen=false),
    a = FitParam(a, lower_limit=0.9, upper_limit=0.998, frozen=false),
    h = FitParam(h, lower_limit=3., upper_limit=100., frozen=false),
    θ = FitParam(θ, lower_limit=1., upper_limit=89., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

LPConv = XS_LampPostJohannsen(;
    E = FitParam(EConv, lower_limit=1., upper_limit=1.05, frozen=false),
    a = FitParam(a, lower_limit=0.9, upper_limit=0.998, frozen=false),
    h = FitParam(h, lower_limit=3., upper_limit=100., frozen=false),
    θ = FitParam(θ, lower_limit=1., upper_limit=89., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(Kₛ, lower_limit = 0., upper_limit=1., frozen=false),
    Γ = FitParam(Γ, lower_limit = 1.2, upper_limit = 1.9, frozen = false),
    A_Fe = FitParam(A_Fe, lower_limit = 0., upper_limit = 5., frozen = true),
    logXi = FitParam(logXi, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(density, lower_limit=1., upper_limit=30., frozen = false), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(ηH, lower_limit=5.0, upper_limit=13.0, frozen=false)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(Kₚ, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=1.2, upper_limit=1.9, frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LP) 

# Defining the model
model = SpectralFitting.Constant(value=FitParam(1.)) * AB * (PL + ConvModel(SP) + LP)

# Defining the problem and binding parameters
prob = FittingProblem(model => data)

# Binding photon indexes
bind!(prob, (1, :a1, :a) => (1, :a2, :Γ))

# Binding line profile parameters
bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination) => (1, :a3, :θ))
bind!(prob, (1, :c1, :a) => ((1, :a3, :a)))
bind!(prob, (1, :c1, :h) => ((1, :a3, :h)))
bind!(prob, (1, :c1, :α13) => ((1, :a3, :α13)))
bind!(prob, (1, :c1, :ϵ3) => ((1, :a3, :ϵ3)))

# Fitting
result = SpectralFitting.fit(
    prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1000)
)

display(result)

time = now()

plt = PlotResult(data, result[1])
display(plt)
savefig(plt, "output/FitResult_$(data.user_data.object)_$time.png")

close("all")

@save "output/FitResult_$(data.user_data.object)_$time.bson" result

open("output/FitResult_$(data.user_data.object)_$time.txt", "w") do io
    JSON3.pretty(io, vars)
end