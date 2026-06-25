using SpectralFitting
using Plots
using LaTeXStrings

include("LampPostModelFit.jl")

default(titlefont = (12, "serif"), 
    guidefont = (12, "serif"), 
    legendfont = (12, "serif"), 
    tickfont = (10, "serif"), 
    minorgrid = true, 
    gridalpha=0.5,
    minorgridalpha=0.3
)

pyplot()

DATADIR = "/home/brad/Documents/SummerInternship/data/"
EXTENSION = "_sr_1000.pha"

function Residuals(result, domain; bounds = (5, 7.5))
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1]
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

function FitPowerLaw(data)

    # Fitting a power law model to the data
    model = PowerLaw()
    prob = FittingProblem(model => data)

    SpectralFitting.fit(prob, LevenbergMarquadt();verbose=true)
end

function FitLineProfile(data, energy, Model=LampPostJohannsen; kwargs...)

    model = Model(E=energy)
    prob = FittingProblem(model => data)

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

# Reading the data
path = joinpath(DATADIR, "$(files[3])B01$(EXTENSION)")
data = LoadData(path)

# Fitting a power law
println("Fitting power law...")
powerLawResult = FitPowerLaw(data)

# Taking the residuals from the power law fit
domain = SpectralFitting.plotting_domain(data)
residuals = Residuals(powerLawResult, domain)

# Plotting
hline([0], c=:grey, label=nothing)
vline!([6.4], c=:grey, linestyle=:dash, label=nothing)

display(
    plot!(residuals, seriestype = :stepmid, c=:black, 
        xlabel="Energy (keV)", ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)", 
        xlims=(5,7.5), legend=:topleft
    )
)

# Fitting a line profile using Gradus' Johannsen metric
println("Fitting Johannsen...")
johannsenResult = FitLineProfile(residuals,
    FitParam(6.4),
    LampPostJohannsen
)

display(plot!(johannsenResult, c=:blue))

# Fitting a line profile using Gradus' Kerr metric
println("Fitting Kerr...")
kerrResult = FitLineProfile(residuals,
    FitParam(6.4),
    LampPostKerr
)

display(plot!(kerrResult, c=:red))
