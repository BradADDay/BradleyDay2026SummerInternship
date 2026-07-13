using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load
using BenchmarkTools

# Loading in files
include("utils/UTILS.jl")

# Setting Filepaths
const ROOT = "/home/brad/Documents/SummerInternship/"
const DATADIR = joinpath(ROOT, "data")
const NUSTAR_EXTENSION = "_sr_1000.pha"
const OUTPUT = joinpath(ROOT, "output/")

# ======================================================================================
# Setup
# ======================================================================================
##
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

dataRange = (3,10)
index = 2

for index in eachindex(files)
    
    # Reading the data
    pathA = joinpath(DATADIR, "$(files[index])A01$(NUSTAR_EXTENSION)")
    dataA = LoadData(pathA; dataRange)
    domainA = SpectralFitting.plotting_domain(dataA)

    pathB = joinpath(DATADIR, "$(files[index])B01$(NUSTAR_EXTENSION)")
    dataB = LoadData(pathB; dataRange)
    domainB = SpectralFitting.plotting_domain(dataB)

    # Allowing the energy to vary between 6.4 keV (neutral/weakly ionised) and 7 (H-like iron)
    energy = FitParam(6.4, lower_limit=6.4, upper_limit=7, frozen=false)

    # ======================================================================================
    # Kerr metric
    # ======================================================================================

    println("Fitting Kerr...")

    # Fitting the table model with the deformation parameters set to 0
    kerrResult = FitPowerLawLineProfile(dataA, dataB; E=copy(energy), α13=FitParam(0.0, frozen=true), ϵ3=FitParam(0.0, frozen=true))

    # Plotting
    PlotFits(dataA, kerrResult[1], dataB, kerrResult[2]; title="Kerr Fit")
    PlotResiduals(dataA, kerrResult[1])

    LP, PL = GetParams(kerrResult; model="K")

    # ======================================================================================
    # Johannsen metric
    # ======================================================================================

    println("Fitting Johannsen...")

    # Fitting the table model
    johannsenResult = FitPowerLawLineProfile(dataA, dataB; E=energy)

    # Plotting
    PlotFits(dataA, johannsenResult[1], dataB, johannsenResult[2]; title="Johannsen Fit")
    PlotResiduals(dataA, johannsenResult[1])

    # ======================================================================================
    # Saving
    # ======================================================================================

    savefig(joinpath(OUTPUT, "TableKerr$(files[index]).png"))
    @save joinpath(OUTPUT, "TableKerrResult$(files[index]).bson") kerrResult
    savefig(joinpath(OUTPUT, "TableJohannsen$(files[index]).png"))
    @save joinpath(OUTPUT, "TableJohannsenResult$(files[index]).bson") johannsenResult

    close("all")
end

CompleteSound()

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