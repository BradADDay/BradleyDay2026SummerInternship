using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load
using BenchmarkTools

# Loading in files
include("FittingModels.jl")
include("ParameterVariations.jl")

# Setting plotting defaults
default(titlefont = (12, "serif"), 
    guidefont = (12, "serif"), 
    legendfont = (12, "serif"), 
    tickfont = (10, "serif"), 
    gridalpha=0.,
    minorticks=true
)
pyplot()

# Setting Filepaths
const ROOT = "/home/brad/Documents/SummerInternship/"
const DATADIR = joinpath(ROOT, "data")
const EXTENSION = "_sr_1000.pha"
const OUTPUT = joinpath(ROOT, "output/")

function CompleteSound()
    """A function to play a sound when the program finishes running"""
    y, fs = wavread("$ROOT/Code/utils/complete.wav")
    wavplay(y, fs)
end

function LoadData(path; dataRange=(3,12))
    """Load in an OGIP dataset from a given path and 
    curtail it to an energy range"""
    # Reading the dataset
    data = OGIPDataset(path)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,dataRange...)
end

function BindParameters(modelA, modelB, dataA, dataB)
    """Bind the parameters of a Johannsen model together between two datasets"""
    prob = FittingProblem(modelA => dataA, modelB => dataB)

    bind!(prob, (1, :a1, :a)   => (2, :a1, :a))
    bind!(prob, (1, :a1, :h)   => (2, :a1, :h))
    bind!(prob, (1, :a1, :θ)   => (2, :a1, :θ))
    bind!(prob, (1, :a1, :E)   => (2, :a1, :E))
    bind!(prob, (1, :a2, :a)   => (2, :a2, :a))
    bind!(prob, (1, :a1, :α13) => (2, :a1, :α13))
    bind!(prob, (1, :a1, :ϵ3)  => (2, :a1, :ϵ3))

    return prob
end

function FitPowerLawLineProfile(dataA, dataB; kwargs...)
    """Fit a composite model of a power law and line profile 
    from the Johannsen table model"""

    # Defining the models for the two datasets
    modelA = XS_LampPostJohannsen(;kwargs...) + PowerLaw()
    modelB = XS_LampPostJohannsen(;kwargs...) + PowerLaw()

    # Binding the parameters together, excluding only the normalisation
    prob = BindParameters(modelA, modelB, dataA, dataB)

    # Fitting the model to the data
    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true)
end

function DualSpectrumPlot(plotA, plotB; bounds=(5,7.5), kwargs...)
    """Plotting the results of the fit to the two datasets"""

    # Plotting
    # Defining the layout
    layout = @layout [a{0.001w} (2,1)]

    # Plotting the y axis label
    yAxis = plot([0], c=:white; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)", 
        framestyle=:none, gridalpha=0, legend=false, xlims=(1,2)
    )

    # Plotting the data
    figure = plot(yAxis, plotA, plotB; layout=layout, link=:x,  
        xlims=bounds, xminorticks=4, margin=1mm, kwargs...
    )

    display(figure)

    return figure
end

function PlotSpectrum(data; xlabel=nothing, ylabel=nothing)
    """Plot a spectrum with a vertical line denoting the iron Kα line"""

    # Plotting the vertical line
    plot = vline([6.4], c=:black, linestyle=:dash, label=nothing)

    # Plotting the spectrum
    plot!(data; seriestype = :stepmid, c=:black, 
        legend=:outerright, framestyle=:box, xminorticks=4,
        xlabel=xlabel, ylabel=ylabel, edgecolor=nothing
    )

    # Functionality to turn off the x ticks
    if isnothing(xlabel)
        plot!(xformatter= _-> "")
    end

    return plot
end

function SeparateModel(model; johannsen=false)
    params = model.u
    K, h, θ, a = params[1:4]
    Kp, ap = params[end-1:end]

    α13, ϵ3 = params[5:6]
    return LampPostJohannsen(;
        K = FitParam(K),
        h = FitParam(h),
        θ = FitParam(θ),
        a = FitParam(a),
        α13 = FitParam(α13),
        ϵ3 = FitParam(ϵ3)
    ), PowerLaw(K = FitParam(Kp), a = FitParam(ap))
end

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
index = 1

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
domainA = SpectralFitting.plotting_domain(dataA)

pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB; dataRange)
domainB = SpectralFitting.plotting_domain(dataB)

# Allowing the energy to vary between 6.4 keV (neutral/weakly ionised) and 7 (H-like iron)
energy = FitParam(6.4, lower_limit=6.4, upper_limit=7)

# ======================================================================================
# Kerr metric
# ======================================================================================

println("Fitting Kerr...")

# Fitting the table model with the deformation parameters set to 0
kerrResult = FitPowerLawLineProfile(dataA, dataB; E=energy, α13=FitParam(0.0, frozen=true), ϵ3=FitParam(0.0, frozen=true))

# Plotting
plotA = PlotSpectrum(dataA)
plot!(plotA, kerrResult[1])

plotB = PlotSpectrum(dataB; xlabel="Energy (keV)")
plot!(plotB, kerrResult[2])

kerrFigure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

# ======================================================================================
# Johannsen metric
# ======================================================================================

println("Fitting Johannsen...")

# Fitting the table model
johannsenResult = FitPowerLawLineProfile(dataA, dataB; E=energy)

# Plotting
plotA = PlotSpectrum(dataA)
plot!(plotA, johannsenResult[1])

plotB = PlotSpectrum(dataB; xlabel="Energy (keV)")
plot!(plotB, johannsenResult[2])

johanFigure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

# ======================================================================================
# Saving
# ======================================================================================

savefig(joinpath(OUTPUT, "TableKerr$(files[index]).png"))
@save joinpath(OUTPUT, "TableKerrResult$(files[index]).bson") kerrResult
savefig(joinpath(OUTPUT, "TableJohannsen$(files[index]).png"))
@save joinpath(OUTPUT, "TableJohannsenResult$(files[index]).bson") johannsenResult

CompleteSound()