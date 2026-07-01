using SpectralFitting
using Plots
using LaTeXStrings
using Measures
using WAV
using BSON: @save, @load

include("FittingModels.jl")
include("ParameterVariations.jl")

default(titlefont = (12, "serif"), 
    guidefont = (12, "serif"), 
    legendfont = (12, "serif"), 
    tickfont = (10, "serif"), 
    gridalpha=0.,
    minorticks=true
)

pyplot()

const ROOT = "/home/brad/Documents/SummerInternship/"
const DATADIR = joinpath(ROOT, "data")
const EXTENSION = "_sr_1000.pha"
const OUTPUT = joinpath(ROOT, "output/")

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

function DefineProblem(modelA, modelB, dataA, dataB)

    prob = FittingProblem(modelA => dataA, modelB => dataB)

    bind!(prob, (1, :a1, :a)   => (2, :a1, :a))
    bind!(prob, (1, :a1, :h)   => (2, :a1, :h))
    bind!(prob, (1, :a1, :θ)   => (2, :a1, :θ))
    bind!(prob, (1, :a1, :E)   => (2, :a1, :E))
    bind!(prob, (1, :a2, :a)   => (2, :a2, :a))

    return prob
end

function FitPowerLawLineProfile(dataA, dataB, energy::FitParam, kerrResult::FitResult; kwargs...)

    LP1 = LampPostJohannsen(K=FitParam(kerrResult[1].u[1]), h=FitParam(kerrResult[1].u[2]), E=energy, θ=FitParam(kerrResult[1].u[4]), a=FitParam(kerrResult[1].u[5]))
    LP2 = LampPostJohannsen(K=FitParam(kerrResult[2].u[1]), h=FitParam(kerrResult[2].u[2]), E=energy, θ=FitParam(kerrResult[2].u[4]), a=FitParam(kerrResult[2].u[5]))
    PL1 = PowerLaw(K=FitParam(kerrResult[1].u[6], frozen=true), a=FitParam(kerrResult[1].u[7], frozen=true))
    PL2 = PowerLaw(K=FitParam(kerrResult[2].u[6], frozen=true), a=FitParam(kerrResult[2].u[7], frozen=true))

    prob = DefineProblem(LP1+PL1, LP2+PL2, dataA, dataB)

    bind!(prob, (1, :a1, :α13) => (2, :a1, :α13))
    bind!(prob, (1, :a1, :ϵ3)  => (2, :a1, :ϵ3))

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, kwargs...)
end

function FitPowerLawLineProfile(dataA, dataB, energy::FitParam; kwargs...)

    modelA = LampPostKerr(E=energy) + PowerLaw()
    modelB = LampPostKerr(E=energy) + PowerLaw()

    prob = DefineProblem(modelA, modelB, dataA, dataB)

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
        xlabel=xlabel, ylabel=ylabel, edgecolor=nothing
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

index = 3
dataRange = (3,10)

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
domainA = SpectralFitting.plotting_domain(dataA)

pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB; dataRange)
domainB = SpectralFitting.plotting_domain(dataB)

## ===============================================================

# Fitting a power law and line profile simultaneously
# allowing the energy to vary between 6.4 keV (neutral/weakly ionised) and 7 (H-like iron)
energy = FitParam(6.4, lower_limit=6.4, upper_limit=7)

pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
domainA = SpectralFitting.plotting_domain(dataA)

pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB; dataRange)

# Kerr metric
println("Fitting Kerr...")
kerrResult = FitPowerLawLineProfile(dataA, dataB, energy)
@save joinpath(OUTPUT, "kerrResult$(files[index]).bson") kerrResult

# Johannsen metric
println("Fitting Johannsen...")
johannsenResult = FitPowerLawLineProfile(dataA, dataB, energy, kerrResult)
@save joinpath(OUTPUT, "johannsenResult$(files[index]).bson") johannsenResult

## ================================================================

@load joinpath(OUTPUT, "kerrResult$(files[index]).bson") kerrResult
@load joinpath(OUTPUT, "johannsenResult$(files[index]).bson") johannsenResult

plotA = PlotSpectrum(dataA)
plot!(plotA, kerrResult[1])

plotB = PlotSpectrum(dataB; xlabel="Energy (keV)")
plot!(plotB, kerrResult[2])

figure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

savefig(joinpath(OUTPUT, "Kerr$(files[index]).png"))

plotA = PlotSpectrum(dataA)
plot!(plotA, johannsenResult[1])

plotB = PlotSpectrum(dataB; xlabel="Energy (keV)")
plot!(plotB, johannsenResult[2])

figure = DualSpectrumPlot(plotA, plotB; bounds=(3,10))

savefig(joinpath(OUTPUT, "Johannsen$(files[index]).png"))