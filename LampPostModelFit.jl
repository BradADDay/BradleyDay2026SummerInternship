using Gradus, SpectralFitting, Plots

# Lamppost model for a Johannsen metric
struct LampPostJohannsen{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation"
    K::T
    "Corona Height"
    h::T
    "Line Energy"
    E::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
    "α13"
    α13::T
    "ϵ3"
    ϵ3::T
end

# Utility function for instantiation
function LampPostJohannsen(;K=1., h=10., E=1., R_in=-1., R_out=400., θ=60.,
                            a=0.998, α13=0., ϵ3=0.)
    K = FitParam(K)
    h = FitParam(h, lower_limit=1.5, upper_limit=50., frozen=false)
    E = FitParam(E, lower_limit=1., upper_limit=10., frozen=false)
    R_in = FitParam(R_in, lower_limit=-Inf, frozen=true)
    R_out = FitParam(R_out, lower_limit=-Inf, frozen=true)
    θ = FitParam(θ, lower_limit=5., upper_limit=85., frozen=false)
    a = FitParam(a, lower_limit=-0.998, upper_limit=0.998, frozen=false)
    α13 = FitParam(α13, lower_limit=0, upper_limit=50., frozen=false)
    ϵ3 = FitParam(ϵ3, lower_limit=0, upper_limit=30., frozen=false)
    LampPostJohannsen(K, h, E, R_in, R_out, θ, a, α13, ϵ3)
end

# Defining how the fit works
function SpectralFitting.invoke!(output, input, model::LampPostJohannsen)

    # Scaling for energy
    domain = copy(input) / model.E

    # Instantiating the metric and observer position
    m = JohannsenMetric(;a = model.a, α13 = model.α13, ϵ3 = model.ϵ3)
    x = SVector(0.0, 1e4, deg2rad(model.θ), 0.0)

    # Setting the inner radius to a multiple of the isco if set below zero
    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    # Instantiating the disk
    d = ThinDisc(R_In, Inf)

    # Generating an emissivity profile
    emissivityModel = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, emissivityModel)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; bins = domain, 
                       method = TransferFunctionMethod(), minrₑ=R_In, 
                       maxrₑ=model.R_out, numrₑ = 30)
    
    output .= flux[1:end-1]
    
end