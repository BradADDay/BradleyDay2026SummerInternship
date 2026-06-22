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

function LoadData(path, plot=true)

    data = OGIPDataset(path)

    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,3.0,75.0)

    if plot
        display(plot(data))
    end

    return data
end

function FitPowerLaw(data)

    model = PowerLaw()
    prob = FittingProblem(model => data)

    result = SpectralFitting.fit(prob, LevenbergMarquadt();verbose=true)

    return result
end

function FitLineProfile(domain, data, bounds, energy)

    # Filtering to include interesting region
    data = data[(domain .> bounds[1]) .& (domain .< bounds[2])]
    domain = domain[(domain .> bounds[1]) .& (domain .< bounds[2])]

    data = InjectiveData(domain, data, name="Data")

    model = LampPostJohannsen(E=energy)
    prob = FittingProblem(model => data)

    SpectralFitting.fit(lineProb, LevenbergMarquadt(); autodiff = :finite, verbose=true)
end

files = ["nu80402315002A01_sr_1000.pha", "nu80402315002B01_sr_1000.pha", "nu80402315004A01_sr_1000.pha", "nu80402315004B01_sr_1000.pha", "nu80402315006A01_sr_1000.pha", "nu80402315006B01_sr_1000.pha", "nu80402315008A01_sr_1000.pha", "nu80402315008B01_sr_1000.pha", "nu80402315010A01_sr_1000.pha", "nu80402315010B01_sr_1000.pha", "nu80402315012A01_sr_1000.pha", "nu80402315012B01_sr_1000.pha", "nu80502304002A01_sr_1000.pha", "nu80502304002B01_sr_1000.pha", "nu80502304004A01_sr_1000.pha", "nu80502304004B01_sr_1000.pha", "nu80502304006A01_sr_1000.pha", "nu80502304006B01_sr_1000.pha"]

path = joinpath(DATADIR, files[1])
data = LoadData(path)

powerLawResult = FitPowerLaw(data)

domain = SpectralFitting.plotting_domain(data)
residuals = Residuals(result)

lineProfileResult = FitLineProfile(residuals, (5.5, 7.5), energy)

plot(domain, residuals)
plot!(lineProfileResult)