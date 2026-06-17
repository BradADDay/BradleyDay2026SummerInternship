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

function LampPostJohannsen(;K = FitParam(1.0),
    h = FitParam(10., lower_limit=1.5, upper_limit=50., frozen=false),
    E = FitParam(1., lower_limit=1., upper_limit=10., frozen=false),
    R_in = FitParam(-1., lower_limit=-Inf, frozen=true),
    R_out = FitParam(100., lower_limit=-Inf, frozen=true),
    θ = FitParam(60., lower_limit=5., upper_limit=85.),
    a = FitParam(0.998, lower_limit=-0.998, upper_limit=1.),
    α13 = FitParam(0., lower_limit=0, upper_limit=50., frozen=false),
    ϵ3 = FitParam(0., lower_limit=0, upper_limit=30., frozen=false))

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

    d = ThinDisc(R_In, model.R_out)

    emissivityModel = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, emissivityModel)

    data = lineprofile(m, x, d, profile; bins = copy(input), method = TransferFunctionMethod(), minrₑ=R_In, maxrₑ=model.R_out, numrₑ = 100)
    output .= data[2][1:end-1]
end

struct LampPost{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Corona Height"
    h::T
    "Line energy"
    E::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
end

function LampPost(;K = FitParam(1.0),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 10., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LampPost(K, h, E, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::LampPost)
    g_domain = copy(domain)
    
    m = KerrMetric(;a = model.a)
    x_obs = SVector(0.0, 1e3, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, model.R_out)

    mode = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, mode)

    data = lineprofile(m, x_obs, d, profile ;bins = g_domain, method = TransferFunctionMethod(), numrₑ = 30)
    output .= data[2][1:end-1]
end