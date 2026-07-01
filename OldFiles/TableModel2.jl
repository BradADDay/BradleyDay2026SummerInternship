using SpectralFitting
using MultiLinearInterpolations
using Interpolations
using Plots

const params = 5

const DATAFILE = "/home/brad/Documents/SummerInternship/Code/model.FITS"

const ROOT = "/home/brad/Documents/SummerInternship/"
const DATADIR = joinpath(ROOT, "data")
const EXTENSION = "_sr_1000.pha"
const OUTPUT = joinpath(ROOT, "output/")

data = TableModelData(Val(params), DATAFILE)
const table = TableModelInterpolation(data)

function LoadData(path; dataRange=(3,12))

    # Reading the dataset
    data = OGIPDataset(path)

    # Regrouping, normalising, dropping bad channels and curtailing
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data,dataRange...)
end

# Lamppost model for a Johannsen metric
struct XS_LampPostJohannsen{T} <: AbstractSpectralModel{T, Additive}
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

# Utility function for instantiation
function XS_LampPostJohannsen(;
    K = FitParam(1.),
    h = FitParam(10., lower_limit=1.5, upper_limit=30., frozen=false),
    E = FitParam(1., lower_limit=1., upper_limit=10., frozen=false),
    θ = FitParam(60., lower_limit=5., upper_limit=85., frozen=false),
    a = FitParam(0.998, lower_limit=-0.998, upper_limit=0.998, frozen=false),
    α13 = FitParam(0., upper_limit=50., frozen=false),
    ϵ3 = FitParam(0., upper_limit=30., frozen=false)
    )

    XS_LampPostJohannsen(K, E, a, h, θ, α13, ϵ3)
end

function SpectralFitting.invoke!(output, input, model::XS_LampPostJohannsen)

    println(model)

    let E = model.E, a = model.a, h = model.h, θ = model.θ, α13 = model.α13, ϵ3 = model.ϵ3
        # Scaling for energy
        domain = copy(input) / E

        open("Code/test.txt", "a") do io
            write(io, "1")
        end

        profile = SpectralFitting.interpolate_table!(table, a, h, θ, α13, ϵ3)
        open("Code/test.txt", "a") do io
            write(io, "2")
        end

        interp = linear_interpolation(table.data.energy_bins[1:end-1], profile, extrapolation_bc=Line())

        open("Code/test.txt", "a") do io
            write(io, "3")
        end

        flux = interp.(domain)

        output .= flux[1:end-1]

        println("TEST")

    end
end

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

index = 1
dataRange = (3,10)

# Reading the data
pathA = joinpath(DATADIR, "$(files[index])A01$(EXTENSION)")
dataA = LoadData(pathA; dataRange)
domainA = SpectralFitting.plotting_domain(dataA)

model = XS_LampPostJohannsen(E=FitParam(6.4))

modelA = model + PowerLaw()

prob = FittingProblem(modelA => dataA)

result = SpectralFitting.fit(prob, LevenbergMarquadt(); autodiff = :finite)

plot(result)
plot!(dataA)
