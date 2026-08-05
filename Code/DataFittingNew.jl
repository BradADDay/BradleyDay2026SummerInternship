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
#= data = LoadData(
    "data/xa201014010xtd_src_rebinned.pi", 
    "data/xa201014010xtd_bgd_rebinned.pi", 
    "data/xa201014010xtd_p0300000a0_src.rmf", 
    "data/xa201014010xtd_p0300000a0_ptsrc.arf";
    dataRange=(4,10)
) =#
data = LoadData(
    "data/xa300049010xtd_src_rebinned.pi", 
    "data/xa300049010xtd_bgd_rebinned.pi", 
    "data/xa300049010xtd_p031300010_src.rmf", 
    "data/xa300049010xtd_p031300010_ptsrc.arf";
    dataRange=(3,10)
)

# Lamp post line profile model
# TODO: Increase range of possible corona heights?

# Line profile
a = 0.998
h = 150.0
θ = 11.514
EConv = 1.0200

# Spectrum
Kₛ = 5.4446e-8
Γ = 3.0939
A_Fe = 2.9024
logXi = 0.044902
density = 18.722

# Photoelectric absorption
ηH = 1.4495

# Power Law
Kₚ = 0.029641

# Gaussian
Kₙ = 8.7610e-5
μ = 7.5489
σ = 0.16087

vars = Dict(
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
    "Kₚ"=>Kₚ,
    "μ"=>μ,
    "σ"=>σ
)

LPConv = XS_LampPostJohannsen(;
    E = FitParam(EConv, lower_limit=1., upper_limit=1.05, frozen=false),
    a = FitParam(a, lower_limit=0.9, upper_limit=0.998, frozen=false),
    h = FitParam(h, lower_limit=3., upper_limit=150., frozen=false),
    θ = FitParam(θ, lower_limit=1., upper_limit=89., frozen=false),
    α13 = FitParam(0., lower_limit=-6., upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-8., upper_limit=10., frozen=true)
) 

# Spectrum model, XILLVER
SP = XillverD5(
    K = FitParam(Kₛ, lower_limit = 0., upper_limit=Inf, frozen=false),
    Γ = FitParam(Γ, lower_limit = 1., upper_limit = 5., frozen = false),
    A_Fe = FitParam(A_Fe, lower_limit = 0., upper_limit = 10., frozen = false),
    logXi = FitParam(logXi, lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(density, lower_limit=1., upper_limit=30., frozen = false), 
    inclination = FitParam(θ, lower_limit=5., upper_limit=85., frozen = false)
)

# Photoelectric absorption model
AB = PhotoelectricAbsorption(
    ηH = FitParam(ηH, lower_limit=0.0, upper_limit=13.0, frozen=false)
)

# Power law model for source emission
PL = PowerLaw(
    K = FitParam(Kₚ, lower_limit=0.0, upper_limit=Inf, frozen=false),
    a = FitParam(Γ, lower_limit=1., upper_limit=5., frozen=false)
)

GS = GaussianLine(
    K = FitParam(Kₙ),
    μ = FitParam(μ, lower_limit=7.1, upper_limit=7.75, frozen=false),
    σ = FitParam(σ, frozen=false)
)

# Setting the line profile as a convolution model
ConvModel = AsConvolution(LPConv) 

# Defining the model
model = AB * (PL + ConvModel(SP) + GS)


# Defining the problem and binding parameters
prob = FittingProblem(model => data)

# Binding photon indexes
bind!(prob, (1, :a1, :a) => (1, :a2, :Γ))
# Binding line profile parameters
bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination))

# Fitting
result = SpectralFitting.fit(
    prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(100)
)

display(result)

time = now()

plt = PlotResult(data, result[1]; ylim=[0,50])
display(plt)
close("all")

# Saving
open("output/FitResult_$(data.user_data.object)_$time.txt", "w") do io
    JSON3.pretty(io, vars)
end
savefig(plt, "output/FitResult_$(data.user_data.object)_$time.png")
@save "output/FitResult_$(data.user_data.object)_$time.bson" result
##

begin
    bins = collect(range(3, 10, 1000))
    flux = 50000 .* invokemodel(bins, model)

    PlotSpectrum(data; xlabel="Energy (keV)", ylabel=L"counts s^{-1} keV$^{-1}$")
    ylims!(0,0.5)
    plot!(bins, flux; c=:red)
end