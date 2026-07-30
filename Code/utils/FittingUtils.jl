using LaTeXStrings
using Measures

# Using some of the utils files for plotting, fit definitions
# This also imports SpectralFitting, WAV, CSV, DataFrames, Gradus, and Interpolations
include("PlottingDefaults.jl")
include("FittingModels.jl")
include("GeneralUtils.jl")

# Setting Filepaths
const ROOT = "/home/brad/Documents/SummerInternship/"
const DATADIR = joinpath(ROOT, "data")
const OUTPUT = joinpath(ROOT, "output/")



"""
    Residuals(result, domain; dataRange=(3,10))

Calculate the residuals of a [`SpectralFitting.FitResult`](@ref).
"""
function Residuals(result, domain; dataRange::Tuple{Real, Real}=(3, 10))

    y = calculate_objective!(result, result.u)
    obj, var = get_objective(result), get_objective_variance(result)
    residuals = @. (obj - y) / sqrt(var)

    # Filtering to include interesting region
    residuals = residuals[(domain .> dataRange[1]) .& (domain .< dataRange[2])]
    domain = domain[(domain .> dataRange[1]) .& (domain .< dataRange[2])]

    # Putting into a data object
    InjectiveData(domain, residuals, name="Residuals")

end

"""
    LoadData(path; dataRange=(3,10))

Load in an OGIP dataset from a given spectrum file that includes the response, ancillary
and background files in its header.
"""
function LoadData(path::String; dataRange::Tuple{Real, Real}=(3,10))

    # Reading the dataset
    data = OGIPDataset(path)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data, dataRange...)

end

"""
    LoadData(spectrum, background, response, ancillary; dataRange=(3,10))

Load an OGIP dataset specifying the spectrum, background, response and ancillary files 
separately.
"""
function LoadData(
    spectrum::String, background::String, response::String, ancillary::String; 
    dataRange::Tuple{Real, Real}=(3,10)
    )

    data = OGIPDataset(
        spectrum; 
        background = background, 
        response = response, 
        ancillary = ancillary
    )

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data, dataRange...)
end

"""
    PlotSpectrum(data; xlabel, ylabel)

Plot a spectrum from a dataset alongside a vertical line at 6.4 keV for the Fe Kα line.

If no `xlabel` is provided then the x ticks are disabled, which is used for subplotting in
[`PlotResiduals`](@ref) and [`DualSpectrumPlot`](@ref).
"""
function PlotSpectrum(
        data::SpectralData; xlabel=nothing, ylabel=nothing, 
        dataRange::Tuple{Real, Real}=(3,10)
    )

    # Plotting the vertical line
    plot = vline([6.4], c=:black, linestyle=:dash, label=nothing, xminorticks=4)

    # Plotting the spectrum
    plot!(
        plot, data; seriestype = :stepmid, c=:black, markercolor=:black, lc=:black, lw=0.5, 
        legend=:topright, framestyle=:box, xlabel=xlabel, ylabel=ylabel, label=nothing, 
        xticks=append!(collect(Float64.(range(dataRange...))), 6.4)
    )

    # Functionality to turn off the x ticks
    if isnothing(xlabel)
        plot!(plot, xformatter= _-> "")
    end

    return plot
end

"""
    PlotSpectrum(data, fit; xlabel, ylabel)

Plot a spectrum from a dataset alongside a fit line
"""
function PlotSpectrum(data::SpectralData, fit; xlabel=nothing, ylabel=nothing)
    plot!(PlotSpectrum(data; xlabel=xlabel, ylabel=ylabel), fit; c=:red, lw=1)
end

"""
    PlotResiduals(data::SpectralData, fit; dataRange=(3, 10))

Plot the residuals between a fit and data.
"""
function PlotResiduals(data::SpectralData, fit; dataRange::Tuple{Real, Real}=(3, 10))
    # Getting and plotting the residuals
    residuals = Residuals(fit, SpectralFitting.plotting_domain(data); dataRange=dataRange)
    resPlot = scatter(
        residuals; markershape=:circle, markersize=1.5, c=:black, msw=0, 
        xlabel="Energy (keV)", label=nothing, ylabel="Residuals", 
        xticks=append!(collect(3.:10.), 6.4)
    )

    # Plotting a horizontal line at 0 and a line for the Fe Kα line
    vline!(resPlot, [6.4]; line=(:black, :dash), label=nothing)
    hline!(resPlot, [0]; c=:black, label=nothing)
end

"""
    PlotResult(data, fit; dataRange)

Plot the result of a fit including the data, fit and residuals.
"""
function PlotResult(data::SpectralData, fit; dataRange::Tuple{Real, Real}=(3,10))
    # Plotting the spectrum and fit
    fitPlot = PlotSpectrum(data, fit; ylabel=L"Flux (counts s$^{-1}$ keV$^{-1}$)")

    resPlot = PlotResiduals(data, fit; dataRange)

    # Combining the two plots
    layout = @layout [a{0.85h}; b{0.15h}]
    plot(fitPlot, resPlot; layout = layout)
end

function FitPowerLawLineProfile(data::SpectralData; maxIter::Int, kwargs...)
    model = XS_LampPostJohannsen(;kwargs...) + PowerLaw()

    prob = FittingProblem(model => data)

    SpectralFitting.fit(
        prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=maxIter
    )
end

function FitPowerLawLineProfile(data::SpectralData, model; maxIter::Int, kwargs...)

    prob = FittingProblem(model => data)

    SpectralFitting.fit(
        prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=maxIter
    )

end

# ==========================================================================================
# NuSTAR data
# ==========================================================================================

function BindParameters(modelA, modelB, dataA::SpectralData, dataB::SpectralData)
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

function FitPowerLawLineProfile(
        dataA::SpectralData, dataB::SpectralData; maxIter=10000::Int, kwargs...
    )
    """Fit a composite model of a power law and line profile 
    from the Johannsen table model"""

    # Defining the models for the two datasets
    modelA = XS_LampPostJohannsen(;kwargs...) + PowerLaw()
    modelB = XS_LampPostJohannsen(;kwargs...) + PowerLaw()

    # Binding the parameters together, excluding only the normalisation
    prob = BindParameters(modelA, modelB, dataA, dataB)

    # Fitting the model to the data
    SpectralFitting.fit(
        prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=maxIter
    )
end

function FitPowerLawLineProfile(
        dataA::SpectralData, dataB::SpectralData, modelA, modelB; maxIter::Int
    )
    """Fit a composite model of a power law and line profile 
    from the Johannsen table model"""

    # Binding the parameters together, excluding only the normalisation
    prob = BindParameters(modelA, modelB, dataA, dataB)

    # Fitting the model to the data
    SpectralFitting.fit(
        prob, LevenbergMarquadt(); autodiff = :finite, verbose=true, maxIter=maxIter
    )
end

function DualSpectrumPlot(
        plotA, plotB; bounds::Tuple{Int64, Int64}=(5,7.5), title="", kwargs...
    )
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

    # display(figure)

    figure
end

function PlotFits(
        dataA::SpectralData, fitA, dataB::SpectralData, fitB; 
        bounds::Tuple{Real, Real}=(3, 10), title=""
    )

    plotA = PlotSpectrum(dataA, fitA)

    plotB = PlotSpectrum(dataB, fitB; xlabel="Energy (keV)")

    DualSpectrumPlot(plotA, plotB; bounds=bounds, title=title)

end

function NUStarFitJohannsen(dataA::SpectralData, dataB::SpectralData, energy::FitParam)
    println("Fitting Johannsen...")

    # Fitting the table model
    johannsenResult = FitPowerLawLineProfile(dataA, dataB; E=copy(energy))

    for (i, data) in enumerate([dataA, dataB])
        PlotResiduals(data, johannsenResult[i])
    end

    # Plotting
    PlotFits(dataA, johannsenResult[1], dataB, johannsenResult[2]; title="Johannsen Fit")

    return johannsenResult
end

function NUStarFitKerr(dataA::SpectralData, dataB::SpectralData, energy::FitParam)

    println("Fitting Kerr...")

    # Fitting the table model with the deformation parameters set to 0
    kerrResult = FitPowerLawLineProfile(
        dataA, dataB; E=copy(energy), α13=FitParam(0.0, frozen=true), 
        ϵ3=FitParam(0.0, frozen=true)
    )
    
    for (i, data) in enumerate([dataA, dataB])
        PlotResiduals(data, kerrResult[i])
    end
    
    # Plotting
    PlotFits(dataA, kerrResult[1], dataB, kerrResult[2]; title="Kerr Fit")

    return kerrResult

end

function FitNUStar(
        file::String; extension::String="_sr_1000.pha", dataRange::Tuple{Real, Real}=(3,10),
        energy::FitParam=FitParam(6.4, lower_limit=6.4, upper_limit=7, frozen=false),
        kerr::Bool=true, johannsen::Bool=true
    )
        
    kerrResult, johannsenResult = nothing, nothing

    # Reading the data
    pathA = joinpath(DATADIR, "$(file)A01$(extension)")
    dataA = LoadData(pathA; dataRange)

    pathB = joinpath(DATADIR, "$(file)B01$(extension)")
    dataB = LoadData(pathB; dataRange)

    if kerr
        kerrResult = NUStarFitKerr(dataA, dataB, energy)
    end

    if johannsen
        johannsenResult = NUStarFitJohannsen(dataA, dataB, energy)
    end

    kerrResult, johannsenResult

end

# ==========================================================================================
# XRISM data
# ==========================================================================================

function FitXRISM(
        spectrum::String, background::String, response::String, ancillary::String; 
        dataRange::Tuple{Int64, Int64}=(3,10), 
        energy::FitParam=FitParam(6.4, lower_limit=6.4, upper_limit=7, frozen=false),
        kerr::Bool=true, johannsen::Bool=true, kwargs...
    )

    kerrResult, johannsenResult = nothing, nothing

    data = LoadData(
        spectrum, 
        background, 
        response, 
        ancillary;
        dataRange=dataRange
    )

    if kerr

        println("Fitting Kerr...")

        kerrResult = FitPowerLawLineProfile(
            data; maxIter=Int(1e3), E=energy, 
            α13=FitParam(0.0, frozen=true), 
            ϵ3=FitParam(0.0, frozen=true), kwargs...)
            
        PlotResiduals(data, kerrResult[1])
    end

    if johannsen

        println("Fitting Johannsen...")

        johannsenResult = FitPowerLawLineProfile(data; maxIter=Int(1e3), E=energy, kwargs...)

        PlotResiduals(data, johannsenResult[1])
    end

    kerrResult, johannsenResult

end

# ==========================================================================================
# Bin
# ==========================================================================================

function GetParams(result; model::String="J")
    """
    Get the fit parameters from a composite fit of a johannsen metric line profile and power law
    """

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

function SeparateModel(result, domain, model::String="J")

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