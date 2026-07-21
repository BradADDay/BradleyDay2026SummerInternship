include("utils/GeneralUtils.jl")
include("utils/DeformUtils.jl")

PlotBounds(-0.998)
scatter!((1,-0.75); c=:red)
using Gradus: isco, JohannsenMetric

m = JohannsenMetric(M=1., a=-0.998, α13=-0.75, ϵ3=1.)

isco(m)

plts = ParameterRegions(10., 10.; a=0.998)
PlotRegion(plts...; c=:viridis)
scatter!((1, -1))