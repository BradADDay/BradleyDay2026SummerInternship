using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load
using BenchmarkTools

# Loading in files
include("utils/FittingUtils.jl")

## =====================================================================================
# NUStar
# ======================================================================================

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

FitNUStar(files[2]; kerr=false)

close("all")

## =====================================================================================
# XRISM
# ======================================================================================

energy = FitParam(6.4, lower_limit=5.5, upper_limit=7, frozen=false)

data = LoadData(
    "data/xa000125000xtd_src.pi", 
    "data/xa000125000xtd_bgd.pi", 
    "data/xa000125000xtd_p031100010_src.rmf", 
    "data/xa000125000xtd_p031100010_ptsrc.arf";
    dataRange=(3,10)
)

LPModel = XS_LampPostJohannsen(;
    K = FitParam(1.),
    E = FitParam(6.4, lower_limit=6., upper_limit=7., frozen=false),
    a = FitParam(0.98, lower_limit=0, upper_limit=0.998, frozen=false),
    h = FitParam(3., lower_limit=3., upper_limit=15., frozen=false),
    θ = FitParam(25., lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., lower_limit=-0.4, upper_limit=10., frozen=true),
    ϵ3 = FitParam(0., lower_limit=-0.4, upper_limit=10., frozen=true)
) 

model = LPModel + PowerLaw()

kerrResult = FitPowerLawLineProfile(
    data, 
    model; 
    maxIter=Int(1e3)
)
    
PlotResiduals(data, kerrResult[1])