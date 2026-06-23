using SpectralFitting
using Plots
using LaTeXStrings
using Measures

include("LampPostModelFit.jl")

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

function Residuals(result, domain; bounds = (5, 7.5))
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    residuals = @. (obj - y) / sqrt(var)

    # Filtering to include interesting region
    residuals = residuals[(domain .> bounds[1]) .& (domain .< bounds[2])]
    domain = domain[(domain .> bounds[1]) .& (domain .< bounds[2])]

    # Putting into a data object
    InjectiveData(domain, residuals, name="Residuals")
end

function LoadData(path)

    # Reading the dataset
    data = OGIPDataset(path)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,3.0,75.0)
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
    bind!(prob, (1, :α13) => (2, :α13))
    bind!(prob, (1, :ϵ3)  => (2, :ϵ3))
    bind!(prob, (1, :E)   => (2, :E))

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, kwargs...)
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

index = 1

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA)
domainA = SpectralFitting.plotting_domain(dataA)
pathB = joinpath(DATADIR, "$(files[index])B01$(EXTENSION)")
dataB = LoadData(pathB)
domainB = SpectralFitting.plotting_domain(dataB)

# Fitting a power law
println("Fitting power law...")
powerLawFit = FitPowerLaw(dataA, dataB)

# Taking the residuals from the power law fit
residualsA = Residuals(powerLawFit[1], domainA)
residualsB = Residuals(powerLawFit[2], domainB)

plotA = hline([0], c=:black, linestyle=:dash, label=nothing)
vline!([6.4], c=:black, linestyle=:dash, label=nothing)
plot!(residualsA; seriestype = :stepmid, c=:black, xformatter= _-> "", 
    legend=:outerright, framestyle=:box, xminorticks=4
)

plotB = hline([0], c=:black, linestyle=:dash, label=nothing)
vline!([6.4], c=:black, linestyle=:dash, label=nothing)
plot!(residualsB; seriestype = :stepmid, c=:black, xlabel="Energy (keV)", 
    legend=:outerright, framestyle=:box, xminorticks=4
)

# Plotting
layout = @layout [a{0.001w} (2,1)]
yAxis = plot([0], c=:white; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)", 
    framestyle=:none, gridalpha=0, legend=false, xlims=(1,2)
)
figure = plot(yAxis, plotA, plotB; layout=layout, link=:x,  
    xlims=(5,7.5), xminorticks=4, margin=0cm
)

display(figure)

# Fitting a line profile using Gradus' Johannsen metric
println("Fitting Johannsen...")
johannsenResult = FitLineProfile(residualsA, residualsB,
    FitParam(6.4),
    LampPostJohannsen; maxIter=5
)

plot!(plotA, johannsenResult[1], c=:blue)
plot!(plotB, johannsenResult[2], c=:blue)

display(figure)

# Fitting a line profile using Gradus' Kerr metric
println("Fitting Kerr...")
kerrResult = FitLineProfile(residuals,
    FitParam(6.4),
    LampPostKerr; maxIter=5
)

plot!(plotA, kerrResult[1], c=:red)
plot!(plotB, kerrResult[2], c=:red)

display(figure)