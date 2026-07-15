using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load
using BenchmarkTools

# Loading in files
include("utils/FittingUtils.jl")

# ======================================================================================
# Setup
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

##

energy = FitParam(6.4, lower_limit=6.4, upper_limit=7, frozen=false)

data = OGIPDataset(
    "data/xa000125000xtd_src.pi"; 
    background = "data/xa000125000xtd_bgd.pi", 
    response = "data/xa000125000xtd_p031100010_src.rmf", 
    ancillary = "data/xa000125000xtd_p031100010_ptsrc.arf"
)

# Regrouping, normalising, dropping bad channels and curtailing
regroup!(data)
normalize!(data)
drop_bad_channels!(data)
mask_energies!(data, 3, 10)

result = FitPowerLawLineProfile(data; E=energy, α13=FitParam(0., frozen=true), ϵ3=FitParam(0., frozen=true))

PlotResiduals(data, result[1])