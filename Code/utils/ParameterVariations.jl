using Gradus: JohannsenMetric, SVector, TransferFunctionMethod, BinningMethod, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, rendergeodesics, ConstPointFunctions
using Colors: diverging_palette
# using Measures
using LaTeXStrings

include("PlottingDefaults.jl")

# =======================================================================
# Functions
# =======================================================================

"""
    RedshiftColormap(img)

Get a colormap for a redshift render using Gradus. Values greater than 1 get coloured 
blue (blueshifting) and less than 1 get assigned red (redshifting). Non-doppler shifted 
regions get coloured white.

See also [`RenderImage`](@ref), [`PlotImage`](@ref).
"""
function RedshiftColormap(img)

    img = collect(Iterators.flatten(img))

    img[isnan.(img)] .= 1

    iMax = maximum(img)
    iMin = minimum(img)

    mid = (1-iMin) / (iMax-iMin)

    cmap = diverging_palette(0, 255, 100, mid=mid, d1=1, d2=1, b=0)

    return cmap, (iMin, iMax)
end

"""
    JohannsenLineProfile(;
        M=1., a=0.998, α13=0., ϵ3=0., α52=0., α22=0., θ=60., h=10., 
        bins = range(0.0, 1.5, 180), minrₑ=-1., maxrₑ=400., numrₑ=100, 
        method=TransferFunctionMethod(), kwargs...
    )

Compute an emission line profile from a Johannsen metric with a lamppost corona and thin 
disk accretion model.

See also [`PlotJohannsenLineProfile`](@ref).

# Arguments
- `method`: The method to be used by Gradus' solver, either `TransferFunctionMethod()` 
(default) or `BinningMethod()`
"""
function JohannsenLineProfile(;
        M=1., a=0.998, α13=0., ϵ3=0., α52=0., α22=0., θ=60., h=10., 
        bins = range(0.0, 1.5, 180), minrₑ=-1., maxrₑ=400., numrₑ=100, 
        method=TransferFunctionMethod(), kwargs...
    )

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

    # Instantiating the metric
    m = JohannsenMetric(M=M, a=a, α13=α13, ϵ3=ϵ3, α52=α52, α22=α22)

    # Setting the inner radius to the ISCO if the entered value is < 0
    if minrₑ < 0.
        minrₑ = isco(m)
    end

    # Disk
    d = ThinDisc(minrₑ, Inf)

    # Pre-computing transfer functions
    tfs = transferfunctions(
        m, x, d;
        maxrₑ= maxrₑ,
        numrₑ=200, 
        minrₑ=minrₑ,
        kwargs...
    )

    # Setting up the model and emissivity profile
    model = LampPostModel(h = h)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    flux = integrate_lineprofile(
        profile, tfs, bins; rmin=minrₑ, rmax=maxrₑ
    )

    flux
end

function RenderImage(;
        M=1., a=0.998, α13=0., ϵ3=0., α52=0., α22=0., θ=60., imageSize=(100,75), 
        αlim=20, βlim=15, kwargs...
    )
    """
    Render a redshift image of the disk
    """

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)
    λ_max = 2x[2]

    # Instantiating the metric
    m = JohannsenMetric(M=M, a=a, α13=α13, ϵ3=ϵ3, α52=α52, α22=α22)

    # Disk
    d = ThinDisc(0.0, 15.0)

    # Redshift point function
    redshift = ConstPointFunctions.redshift(m, x)
    redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()

    # Rendering the image
    α, β, image = rendergeodesics(
        m, x, d, λ_max, pf = redshiftGeometry,
        # image parameters
        image_width = imageSize[1], image_height = imageSize[2],
        αlims = (-αlim, αlim), βlims = (-βlim, βlim), verbose = true
    )

    α, β, image
end

function PlotImage(;
        M=1., a=0.998, α13=0., ϵ3=0., α52=0., α22=0., θ=60., imageSize=(100,75), 
        αlim=20, βlim=15, title="", kwargs...
    )

    α, β, image = RenderImage(; 
        M, a, α13, ϵ3, α52, α22, θ, imageSize, βlim, αlim, kwargs...
    )

    # Generating a colormap accurately depicting red/blueshift
    # cmap, clims = RedshiftColormap(image)
    clims = (0.024515552807013563, 1.2349754807819402)
    
    mid = (1-clims[1]) / (clims[2]-clims[1])

    cmap = diverging_palette(0, 255, 100, mid=mid, d1=1, d2=1, b=0)

    # Looping through each config and plotting/storing the plot
    hmp = heatmap(
        α, β, image, aspect_ratio = 1; title=title, 
        cmap=cmap, clims=clims, cbartitle="Redshift", 
        xlims=(-αlim, αlim), ylims=(-βlim, βlim),
        xlabel=L"\alpha", ylabel=L"\beta"
    )

    α, β, image, hmp

end

function PlotJohannsenLineProfile(;
        M=1., a=0.998, α13=0., ϵ3=0., α52=0., α22=0., θ=60., h=10., E=1.,
        bins=range(0.0, 1.5, 180), minrₑ=-1., maxrₑ=400., numrₑ=100, 
        method=TransferFunctionMethod(), title="", kwargs...
    )

    # Computing the line profile
    flux = JohannsenLineProfile(;
        M, a, α13, ϵ3, α52, α22, θ, h, bins, minrₑ, maxrₑ, numrₑ, method, kwargs...
    )

    # Plotting the line profile
    plt = plot(bins*E, flux; xlabel="Energy (keV)", ylabel="Flux (Arbitrary)", title=title)

    flux, plt
 
end