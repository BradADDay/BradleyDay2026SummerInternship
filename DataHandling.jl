using FITSIO
using Plots

include("LampPostModelFit.jl")
include("ParameterVariations.jl")
include("Defaults.jl")

file = FITS("nu80502304002B01_sr_1000.pha")

data = read(file[2], "COUNTS")
channel = read(file[2], "CHANNEL")

plot(channel, data)