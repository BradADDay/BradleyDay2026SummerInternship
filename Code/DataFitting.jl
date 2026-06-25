using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load

include("LampPostModelFit.jl")
include("ParameterVariations.jl")

default(titlefont = (12, "serif"), 
    guidefont = (12, "serif"), 
    legendfont = (12, "serif"), 
    tickfont = (10, "serif"), 
    gridalpha=0.,
    minorticks=true
)

pyplot()

DATADIR = "/home/brad/Documents/SummerInternship/data/"
EXTENSION = "_sr_1000.pha"
IMAGEDIR = "/home/brad/Documents/SummerInternship/Images/"

function CompleteSound()
    y, fs = wavread(raw"./complete.wav")
    wavplay(y, fs)
end

function Residuals(result, domain, scope; bounds = (5, 7.5))
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    residuals = @. (obj - y) / sqrt(var)

    # Filtering to include interesting region
    residuals = residuals[(domain .> bounds[1]) .& (domain .< bounds[2])]
    domain = domain[(domain .> bounds[1]) .& (domain .< bounds[2])]

    # Putting into a data object
    InjectiveData(domain, residuals, name="Residuals$scope")
end

function LoadData(path; dataRange=(3,12))

    # Reading the dataset
    data = OGIPDataset(path)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,dataRange...)
end

function FitPowerLaw(dataA, dataB)

    # Fitting a power law model to the data
    modelA = PowerLaw()
    modelB = PowerLaw()

    prob = FittingProblem(modelA => dataA, modelB => dataB)
    bind!(prob, (1, :a) => (2, :a))

    SpectralFitting.fit(prob, LevenbergMarquadt();verbose=true)
end

function FitLineProfile(dataA, dataB, energy, Model=LampPostJohannsen; kwargs...)

    modelA = Model(E=energy)
    modelB = Model(E=energy)

    prob = FittingProblem(modelA => dataA, modelB => dataB)
    bind!(prob, (1, :a)   => (2, :a))
    bind!(prob, (1, :h)   => (2, :h))
    bind!(prob, (1, :θ)   => (2, :θ))
    bind!(prob, (1, :E)   => (2, :E))
    if Model == LampPostJohannsen
        bind!(prob, (1, :α13) => (2, :α13))
        bind!(prob, (1, :ϵ3)  => (2, :ϵ3))
    end

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, kwargs...)
end

function FitPowerLawLineProfile(dataA, dataB, energy, Model=LampPostJohannsen; kwargs...)

    modelA = Model(E=energy) + PowerLaw()
    modelB = Model(E=energy) + PowerLaw()

    prob = FittingProblem(modelA => dataA, modelB => dataB)
    bind!(prob, (1, :a1, :a)   => (2, :a1, :a))
    bind!(prob, (1, :a1, :h)   => (2, :a1, :h))
    bind!(prob, (1, :a1, :θ)   => (2, :a1, :θ))
    bind!(prob, (1, :a1, :E)   => (2, :a1, :E))
    bind!(prob, (1, :a2, :a)   => (2, :a2, :a))
    if Model == LampPostJohannsen
        bind!(prob, (1, :a1, :α13) => (2, :a1, :α13))
        bind!(prob, (1, :a1, :ϵ3)  => (2, :a1, :ϵ3))
    end

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, kwargs...)
end

function DualSpectrumPlot(plotA, plotB; bounds=(5,7.5), kwargs...)
    # Plotting
    layout = @layout [a{0.001w} (2,1)]
    yAxis = plot([0], c=:white; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)", 
        framestyle=:none, gridalpha=0, legend=false, xlims=(1,2)
    )
    figure = plot(yAxis, plotA, plotB; layout=layout, link=:x,  
        xlims=bounds, xminorticks=4, margin=1mm, kwargs...
    )

    display(figure)

    return figure
end

function PlotSpectrum(data; xlabel=nothing, ylabel=nothing)
    plot = hline([0], c=:black, linestyle=:dash, label=nothing)
    vline!([6.4], c=:black, linestyle=:dash, label=nothing)
    plot!(data; seriestype = :stepmid, c=:black, 
        legend=:outerright, framestyle=:box, xminorticks=4,
        xlabel=xlabel, ylabel=ylabel
    )
    if isnothing(xlabel)
        plot!(xformatter= _-> "")
    end

    return plot
end

function SeparateModel(model; johannsen=false)
    params = model.u
    K, h, θ, a = params[1:4]
    Kp, ap = params[end-1:end]

    if johannsen
        α13, ϵ3 = params[5:6]
        return LampPostJohannsen(;
            K = FitParam(K),
            h = FitParam(h),
            θ = FitParam(θ),
            a = FitParam(a),
            α13 = FitParam(α13),
            ϵ3 = FitParam(ϵ3)
        ), PowerLaw(K = FitParam(Kp), a = FitParam(ap))
    else
        return LampPostKerr(;
            K = FitParam(K),
            h = FitParam(h),
            θ = FitParam(θ),
            a = FitParam(a)
        ), PowerLaw(K = FitParam(Kp), a = FitParam(ap))
    end
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

index = 2
dataRange = (3,10)

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
@save "bsonFiles/dataA.bson" dataA
domainA = SpectralFitting.plotting_domain(dataA)

pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB; dataRange)
@save "bsonFiles/dataB.bson" dataB
domainB = SpectralFitting.plotting_domain(dataB)

# ================================================================

# Fitting a power law and line profile simultaneously

# allowing the energy to vary between 6.4 keV (neutral/weakly ionised) and 7 (H-like iron)
energy = FitParam(6.4, lower_limit=6.4, upper_limit=7)

# Johannsen metric
println("Fitting Johannsen...")
johannsenResult = FitPowerLawLineProfile(dataA, dataB, energy, LampPostJohannsen)
@save "bsonFiles/johannsenResult$(files[index]).bson" johannsenResult

# Kerr metric
println("Fitting Kerr...")
kerrResult = FitPowerLawLineProfile(dataA, dataB, energy, LampPostKerr)
@save "bsonFiles/kerrResult$(files[index]).bson" kerrResult

# =================================================================
##
@load "bsonFiles/kerrResultnu80402315004.bson" kerrResult
@load "bsonFiles/johannsenResultnu80402315004.bson" johannsenResult

plotA = PlotSpectrum(dataA)
plot!(plotA, kerrResult[1])

plotB = PlotSpectrum(dataB)
plot!(plotB, kerrResult[2]; xlabel="Energy (keV)")

figure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

savefig("$(IMAGEDIR)Kerr$(files[index]).png")

plotA = PlotSpectrum(dataA)
plot!(plotA, johannsenResult[1])

plotB = PlotSpectrum(dataB)
plot!(plotB, johannsenResult[2]; xlabel="Energy (keV)")

figure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

savefig("$(IMAGEDIR)Johannsen$(files[index]).png")