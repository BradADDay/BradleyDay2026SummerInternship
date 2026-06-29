begin
    using Downloads: download
    using FITSFiles
    using DataFramesMeta: DataFrame, @rsubset
    using SpectralFitting
    using BenchmarkTools
end

include("ParameterVariations.jl")
include("Defaults.jl")

"""
h   ∈ [1.5, 30]
θ   ∈ [5, 85]
a   ∈ [-0.998, 0.998]
α13 ∈ [-(1+sqrt(1 - a^2))^3, 50]
ϵ3  ∈ [-(1+sqrt(1 - a^2))^3, 30]

approx 1.4179 per parameter combo
39 hours for 10 iterations of each parameter
"""

num = 2

hs   = range( 1.5  , 30.   , num)
as   = range(-0.998, 00.998, num)
θs   = range( 5.   , 85.   , num)
α13s = range(-1.   , 50.   , num)
ϵ3s  = range(-1.   , 50.   , num)

Card("MODLNAME", "JMetricLP", "The name of the model")

PRIMARY = HDU(
    Primary,
    missing,
    [
        Card("MODLNAME", "JMetricLP", "The name of the model"),
        Card("MODLUNIT", "counts", "The units for the model"),
        Card("REDSHIFT", false, "Whether redshift is to be a parameter"),
        Card("ESCALE", true, "Whether escale is to be a parameter"),
        Card("ADDMODEL", true, "Whether this is an additive model"),
        Card("LOELIMIT", 0, "The model value for energies below those tabulated"),
        Card("HIELIMIT", 0, "The model value for energies above those tabulated"),
        Card("HDUCLASS", "OGIP"),
        Card("HDUCLAS1", "XSPEC TABLE MODEL"),
        Card("HDUVERS", "1.1.0")
    ]
)

PARAMETERS = HDU(
    Bintable,
    [
        Card("TTYPE1", "NAME"),
        Card("TFORM1", "12A"),
        Card("TTYPE2", "METHOD"),
        Card("TFORM2", "J")
    ]
)

data = FITSFiles.fits("/home/brad/Documents/SummerInternship/data/kerrtable.fits")

