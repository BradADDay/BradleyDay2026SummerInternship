using SpectralFitting
using MultiLinearInterpolations
using Interpolations
using Plots
using Gradus

# ===============================================================
# Lamp post table model
# ===============================================================

# The number of parameters the model has
const params = 5

# The path to the model
MODELDATAFILE = "/home/brad/Documents/SummerInternship/Code/utils/model.FITS"

# Reading in the table model and setting it up as an interpolation object
data = TableModelData(Val(params), MODELDATAFILE)
const table = TableModelInterpolation(data)

# Lamppost table model for a Johannsen metric
struct XS_LampPostJohannsen{D, T} <: AbstractTableModel{T, Additive}
    table::D
    "Normalisation"
    K::T
    "Line Energy"
    E::T
    "Spin"
    a::T
    "Corona Height"
    h::T
    "Inclination"
    θ::T
    "α13"
    α13::T
    "ϵ3"
    ϵ3::T
end

# Default constructor
function XS_LampPostJohannsen(;
    K = FitParam(1.),
    E = FitParam(1., lower_limit=1., upper_limit=10., frozen=false),
    a = FitParam(0.998, lower_limit=-0.998, upper_limit=0.998, frozen=false),
    h = FitParam(10., lower_limit=1.5, upper_limit=30., frozen=false),
    θ = FitParam(60., lower_limit=5., upper_limit=85., frozen=false),
    α13 = FitParam(0., upper_limit=50., frozen=false),
    ϵ3 = FitParam(0., upper_limit=30., frozen=false)
    )

    XS_LampPostJohannsen(table, K, E, a, h, θ, α13, ϵ3)
end

# Fitting process
function SpectralFitting.invoke!(output, input, model::XS_LampPostJohannsen)

    # Pulling the necessary parameters
    let table = model.table, E = model.E, a = model.a, h = model.h, θ = model.θ, α13 = model.α13, ϵ3 = model.ϵ3
        
        # Scaling for energy
        domain = copy(input) / E

        # Getting the line profile shape from the table
        profile = SpectralFitting.interpolate_table!(table, a, h, θ, α13, ϵ3)

        # Interpolating the profile to get the result for the input domain
        interp = linear_interpolation(table.data.energy_bins[1:end-1], profile, extrapolation_bc=Line())
        flux = interp.(domain)

        output .= flux[1:end-1]
    end
end

function SpectralFitting.invokemodel!(input, model::XS_LampPostJohannsen)

    output = zeros(length(input)-1)

    SpectralFitting.invoke!(output, input, model)

end

# ===============================================================
# Johannsen Metric
# ===============================================================

function DeformationConstraints(a)
    -(1+sqrt(1 - a^2))^3
end

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
function LampPostJohannsen(;
    K = FitParam(1.),
    h = FitParam(10., lower_limit=1.5, upper_limit=30., frozen=false),
    E = FitParam(1., lower_limit=1., upper_limit=10., frozen=false),
    R_in = FitParam(-1., lower_limit=-Inf, frozen=true),
    R_out = FitParam(400., lower_limit=-Inf, frozen=true),
    θ = FitParam(60., lower_limit=5., upper_limit=85., frozen=false),
    a = FitParam(0.998, lower_limit=-0.998, upper_limit=0.998, frozen=false),
    α13 = FitParam(-1., upper_limit=50., frozen=false),
    ϵ3 = FitParam(-1., upper_limit=30., frozen=false)
    )

    LampPostJohannsen(K, h, E, R_in, R_out, θ, a, α13, ϵ3)
end

# Defining how the fit works
function SpectralFitting.invoke!(output, input, model::LampPostJohannsen)

    # Scaling for energy
    domain = copy(input) / model.E

    # Instantiating the metric, disk and observer position
    m = JohannsenMetric(;a = model.a, α13 = model.α13, ϵ3 = model.ϵ3)
    x = SVector(0.0, 1e4, deg2rad(model.θ), 0.0)
    d = ThinDisc(0., Inf)

    # Setting the inner radius to a multiple of the isco if set below zero
    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    if model.α13 < DeformationConstraints(model.a)
        model.α13 = DeformationConstraints(model.a)
    end

    if model.ϵ3 < DeformationConstraints(model.a)
        model.ϵ3 = DeformationConstraints(model.a)
    end

    # Generating an emissivity profile
    emissivityModel = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, emissivityModel)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; bins = domain, 
                       method = TransferFunctionMethod(), minrₑ=R_In, 
                       maxrₑ=model.R_out, numrₑ = 30)
    
    output .= flux[1:end-1]
    
end

# ===============================================================
# Kerr Metric
# ===============================================================

# Lamppost model for a Johannsen metric
struct LampPostKerr{T} <: AbstractSpectralModel{T, Additive}
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
end

# Utility function for instantiation
function LampPostKerr(;
    K = FitParam(1.),
    h = FitParam(10., lower_limit=1.5, upper_limit=30., frozen=false),
    E = FitParam(1., lower_limit=1., upper_limit=10., frozen=false),
    R_in = FitParam(-1., lower_limit=-Inf, frozen=true),
    R_out = FitParam(400., lower_limit=-Inf, frozen=true),
    θ = FitParam(60., lower_limit=5., upper_limit=85., frozen=false),
    a = FitParam(0.998, lower_limit=-0.998, upper_limit=0.998, frozen=false)
    )

    LampPostKerr(K, h, E, R_in, R_out, θ, a)
end

# Defining how the fit works
function SpectralFitting.invoke!(output, input, model::LampPostKerr)

    # Scaling for energy
    domain = copy(input) / model.E

    # Instantiating the metric, disk and observer position
    m = KerrMetric(;a = model.a)
    x = SVector(0.0, 1e4, deg2rad(model.θ), 0.0)
    d = ThinDisc(0., Inf)

    # Setting the inner radius to a multiple of the isco if set below zero
    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    # Generating an emissivity profile
    emissivityModel = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, emissivityModel)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; bins = domain, 
                       method = TransferFunctionMethod(), minrₑ=R_In, 
                       maxrₑ=model.R_out, numrₑ = 30)
    
    output .= flux[1:end-1]
    
end