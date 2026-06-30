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

data = FITSFiles.fits("test.FITS")

