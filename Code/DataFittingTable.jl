using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load
using BenchmarkTools
using PairPlots

pyplot()

cm2in(x) = 0.3937008x
cm2px(x) = Int(round(cm2in(x) * 100))

# Loading in files
include("FittingModels.jl")
include("ParameterVariations.jl")

# Setting plotting defaults
#= default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (5, "serif"), 
    tickfont = (8, "serif"), 
    gridalpha=0.,
    minorticks=true,
    dpi=300,
    size=cm2px.((8,5))
) =#

default(titlefont = (12, "serif"), 
    guidefont = (10, "serif"), 
    legendfont = (8, "serif"), 
    tickfont = (8, "serif"), 
    gridalpha=0.,
    minorticks=true,
    dpi=300,
    size=cm2px.((12,9))
)

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

function DualSpectrumPlot(plotA, plotB; bounds=(5,7.5), title="", kwargs...)
    """Plotting the results of the fit to the two datasets"""

    # Plotting
    # Defining the layout
    layout = @layout [a{0.001h}; [b{0.001w} (2,1)]]

    Title = plot([0], c=:white; title=title, 
        framestyle=:none, gridalpha=0, legend=false, xlims=(1,2)
    )

    # Plotting the y axis label
    yAxis = plot([0], c=:white; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)", 
        framestyle=:none, gridalpha=0, legend=false, xlims=(1,2)
    )

    # Plotting the data
    figure = plot(Title, yAxis, plotA, plotB; layout=layout, link=:x,  
        xlims=bounds, xminorticks=4, margin=1mm, kwargs...
    )

    display(figure)

    return figure
end

function PlotSpectrum(data::SpectralData; xlabel=nothing, ylabel=nothing)
    """Plot a spectrum with a vertical line denoting the iron Kα line"""

    # Plotting the vertical line
    plot = vline([6.4], c=:black, linestyle=:dash, label=nothing, xminorticks=4)

    # Plotting the spectrum
    plot!(plot, data; seriestype = :stepmid, c=:black, 
        legend=:topright, framestyle=:box,
        xlabel=xlabel, ylabel=ylabel, label=nothing, 
        markercolor=:black, xticks=append!(collect(3.:10.), 6.4),
        lc=:black, lw=0.5
    )

    # Functionality to turn off the x ticks
    if isnothing(xlabel)
        plot!(plot, xformatter= _-> "")
    end

    return plot
end

function PlotSpectrum(data::SpectralData, fit; xlabel=nothing, ylabel=nothing)
    """Plot a spectrum with a vertical line denoting the iron Kα line"""

    plot!(PlotSpectrum(data; xlabel=xlabel, ylabel=ylabel), fit; c=:red, lw=1)
end

function PlotFits(dataA, fitA, dataB, fitB; bounds=(3, 10), title="")

    plotA = PlotSpectrum(dataA, fitA)

    plotB = PlotSpectrum(dataB, fitB; xlabel="Energy (keV)")

    DualSpectrumPlot(plotA, plotB; bounds=bounds, title=title)

end

function SeparateModel(result, domain, model="J")

    values = result.u

    K, E, a, h, θ = values[1:5]
    K2, a2 = values[end-1:end]

    if model == "J"
        α13, ϵ3 = values[6:7]
        LP = XS_LampPostJohannsen(;K=K, E=E, a=a, h=h, θ=θ, α13=α13, ϵ3=ϵ3)
    else
        LP = XS_LampPostJohannsen(;K=K, E=E, a=a, h=h, θ=θ)
    end
    
    println(K2)

    PL = invokemodel(domain, PowerLaw(K = FitParam(K2), a=FitParam(a2)))

    return LP, PL

end

function FitContour(result, params, param1, param2)

    values = copy(result.u)

    stats = zeros(length(param1), length(param2))

    for i in eachindex(param1)
        for j in eachindex(param2)
            values[params[1]] = param1[i]
            values[params[2]] = param2[j]
            stats[i,j] = measure(ChiSquared(), result, values)
        end
    end

    pairplot(stats)
    
    scatter!([result.u[params[1]]], [result.u[params[2]]])
end

#= function FitContour(result, params, param1, param2)

    values = copy(result.u)

    stats = zeros(length(param1), length(param2))

    for i in eachindex(param1)
        for j in eachindex(param2)
            values[params[1]] = param1[i]
            values[params[2]] = param2[j]
            stats[i,j] = measure(ChiSquared(), result, values)
        end
    end

    println(size(stats))

    # 1, 2, and 3 sigma contours
    stdev = std(stats)
    contour(
        param1,
        param2,
        stats .- sum(result.stats),
        levels = [1stdev, 2stdev, 3stdev],
        xlabel = params[3],
        ylabel = params[4]
    )
    scatter!([result.u[params[1]]], [result.u[params[2]]])
end =#

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

dataRange = (3,10)
index = 2

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
domainA = SpectralFitting.plotting_domain(dataA)

pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB; dataRange)
domainB = SpectralFitting.plotting_domain(dataB)

# Allowing the energy to vary between 6.4 keV (neutral/weakly ionised) and 7 (H-like iron)
energy = FitParam(6.4, lower_limit=6.4, upper_limit=7, frozen=false)

# ======================================================================================
# Kerr metric
# ======================================================================================

println("Fitting Kerr...")

# Fitting the table model with the deformation parameters set to 0
kerrResult = FitPowerLawLineProfile(dataA, dataB; E=energy, α13=FitParam(0.0, frozen=true), ϵ3=FitParam(0.0, frozen=true))

# Plotting
PlotFits(dataA, kerrResult[1], dataB, kerrResult[2]; title="Kerr Fit")

FitContour(kerrResult[2], (3, 4, "a", "h"), range(0, 0.998, 50), range(0, 30, 50))

## ======================================================================================
# Johannsen metric
# ======================================================================================

println("Fitting Johannsen...")

# Fitting the table model
johannsenResult = FitPowerLawLineProfile(dataA, dataB; E=energy)

# Plotting
PlotFits(dataA, johannsenResult[1], dataB, johannsenResult[2]; title="Johannsen Fit")

FitContour(johannsenResult[2], (2, 4, "E", "h"), range(1., 10., 50), range(0, 35, 50))

# ======================================================================================
# Saving
# ======================================================================================

savefig(joinpath(OUTPUT, "TableKerr$(files[index]).png"))
@save joinpath(OUTPUT, "TableKerrResult$(files[index]).bson") kerrResult
savefig(joinpath(OUTPUT, "TableJohannsen$(files[index]).png"))
@save joinpath(OUTPUT, "TableJohannsenResult$(files[index]).bson") johannsenResult

CompleteSound()

