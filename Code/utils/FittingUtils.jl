using SpectralFitting
using LaTeXStrings
using Measures

include("PlottingDefaults.jl")
include("FittingModels.jl")
include("GeneralUtils.jl")

function Residuals(result, domain; bounds = (3, 10))
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

function PlotSpectrum(data; xlabel=nothing, ylabel=nothing)
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

function GetParams(result; model="J")

    values = result.u

    K, E, a, h, θ = values[1:5]
    K2, a2 = values[end-1:end]

    LP = (K=K, E, a, h, θ)
    PL = (K=K2, a=a2)

    if model == "J"
        α13, ϵ3 = values[6:7]
        setindex(LP, "α13", α13)
        setindex(LP, "ϵ3", ϵ3)
    end

    LP, PL
end

function SeparateModel(result, domain, model="J")

    LPParams, PLParams = GetParams(result; model=model)

    LP = XS_LampPostJohannsen(;LPParams...)

    PL = PowerLaw(PLParams...)

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
end

function PlotResiduals(data, fit)

    fitPlot = PlotSpectrum(data, fit; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)")

    residuals = Residuals(fit, SpectralFitting.plotting_domain(data))
    resPlot = scatter(residuals; markershape=:circle, markersize=2, c=:black, msw=0, xlabel="Energy (keV)", label=nothing, ylabel="Residuals")
    vline!(resPlot, [6.4]; line=(:black, :dash), label=nothing)
    hline!(resPlot, [0]; c=:black, label=nothing)

    layout = @layout [a{0.75h}; b{0.25h}]

    figure = plot(fitPlot, resPlot; layout = layout)

    display(figure)
    
end

function FitPowerLawLineProfile(data; kwargs...)
    model = XS_LampPostJohannsen(;kwargs...) + PowerLaw()

    prob = FittingProblem(model => data)

    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1e4))
end

# =================================================================================
# NuSTAR data
# =================================================================================

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
    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1e4))
end

function FitPowerLawLineProfile(dataA, dataB, modelA, modelB)
    """Fit a composite model of a power law and line profile 
    from the Johannsen table model"""

    # Binding the parameters together, excluding only the normalisation
    prob = BindParameters(modelA, modelB, dataA, dataB)

    # Fitting the model to the data
    SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=Int(1e4))
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

function PlotFits(dataA, fitA, dataB, fitB; bounds=(3, 10), title="")

    plotA = PlotSpectrum(dataA, fitA)

    plotB = PlotSpectrum(dataB, fitB; xlabel="Energy (keV)")

    DualSpectrumPlot(plotA, plotB; bounds=bounds, title=title)

end