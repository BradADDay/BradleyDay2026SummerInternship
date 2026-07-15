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

kerrResult, johannsenResult = FitXRISM(
    "data/xa000125000xtd_src.pi", 
    "data/xa000125000xtd_bgd.pi", 
    "data/xa000125000xtd_p031100010_src.rmf", 
    "data/xa000125000xtd_p031100010_ptsrc.arf";
    johannsen=false,
    K=FitParam(10.)
)