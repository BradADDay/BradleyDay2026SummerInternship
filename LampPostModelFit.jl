using Gradus, SpectralFitting, Plots

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

function LampPostJohannsen(; K=1., h=10., E=1., R_in=-1., R_out=400., θ=60.,
                            a=0.998, α13=0., ϵ3=0.)
    K = FitParam(K)
    h = FitParam(h, lower_limit=1.5, upper_limit=50., frozen=false)
    E = FitParam(E, lower_limit=1., upper_limit=10., frozen=false)
    R_in = FitParam(R_in, lower_limit=-Inf, frozen=true)
    R_out = FitParam(R_out, lower_limit=-Inf, frozen=true)
    θ = FitParam(θ, lower_limit=5., upper_limit=85.)
    a = FitParam(a, lower_limit=-0.998, upper_limit=0.998)
    α13 = FitParam(α13, lower_limit=0, upper_limit=50., frozen=false)
    ϵ3 = FitParam(ϵ3, lower_limit=0, upper_limit=30., frozen=false)
    LampPostJohannsen(K, h, E, R_in, R_out, θ, a, α13, ϵ3)
end

function SpectralFitting.invoke!(output, input, model::LampPostJohannsen)

    m = JohannsenMetric(;a = model.a, α13 = model.α13, ϵ3 = model.ϵ3)
    x = SVector(0.0, 1e4, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, Inf)

    emissivityModel = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, emissivityModel)

    data = lineprofile(m, x, d, profile; bins = copy(input), 
                       method = TransferFunctionMethod(), minrₑ=R_In, 
                       maxrₑ=model.R_out, numrₑ = 30)

    output .= data[2][1:end-1]
    
end