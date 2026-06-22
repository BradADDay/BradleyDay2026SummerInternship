using SpectralFitting
using Plots

include("LampPostModelFit.jl")

plotly()

DATADIR = "/home/brad/Documents/SummerInternship/data/"

function Residuals(result)
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1]
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    @. (obj - y) / sqrt(var)
end

function LoadData(path)

    data = OGIPDataset(path)

    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,3.0,75.0)

    return data
end

function FitPowerLaw(data)

    model = PowerLaw()
    prob = FittingProblem(model => data)

    SpectralFitting.fit(prob, LevenbergMarquadt();verbose=true)
end

function FitLineProfile(data, energy)

    model = LampPostJohannsen(E=energy)
    prob = FittingProblem(model => data)

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=50)
end

# List of available datasets
files = ["nu80402315002A01_sr_1000.pha", "nu80402315002B01_sr_1000.pha", "nu80402315004A01_sr_1000.pha", "nu80402315004B01_sr_1000.pha", "nu80402315006A01_sr_1000.pha", "nu80402315006B01_sr_1000.pha", "nu80402315008A01_sr_1000.pha", "nu80402315008B01_sr_1000.pha", "nu80402315010A01_sr_1000.pha", "nu80402315010B01_sr_1000.pha", "nu80402315012A01_sr_1000.pha", "nu80402315012B01_sr_1000.pha", "nu80502304002A01_sr_1000.pha", "nu80502304002B01_sr_1000.pha", "nu80502304004A01_sr_1000.pha", "nu80502304004B01_sr_1000.pha", "nu80502304006A01_sr_1000.pha", "nu80502304006B01_sr_1000.pha"]

# Reading the data
path = joinpath(DATADIR, files[2])
data = LoadData(path)

# Fitting a power law
println("Fitting power law...")
powerLawResult = FitPowerLaw(data)

# Taking the residuals from the power law fit
domain = SpectralFitting.plotting_domain(data)
residuals = Residuals(powerLawResult)

# Filtering to include interesting region
bounds = (5, 7.5)
residuals = data[(domain .> bounds[1]) .& (domain .< bounds[2])]
domain = domain[(domain .> bounds[1]) .& (domain .< bounds[2])]

data = InjectiveData(domain, residuals, name="Data")

# Fitting a line profile using Gradus
println("Fitting line profile...")
lineProfileResult = FitLineProfile(data, 6.4)

# Plotting
plot(domain, residuals, seriestype = :stepmid)
plot!(lineResult)