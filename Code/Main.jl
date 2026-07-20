include("utils/PlottingDefaults.jl")

using Gradus

bins, a, h, θ, α13, ϵ3 = (range(0, 3, 200), -0.998, 3.0, 45.0, -0.99, 1.)

println("pos")
x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

println("metric")
# Instantiating the metric
m = JohannsenMetric(M=1., a=a, α13=α13, ϵ3=ϵ3)

println("disk")
minrₑ = Gradus.isco(m)

d = ThinDisc(minrₑ, Inf)

println("lpmodel")
# Setting up the model and emissivity profile

if h < minrₑ
    h = minrₑ
end

model = LampPostModel(h = 10.)
profile = emissivity_profile(m, d, model)

println("lp")
# Computing the line profile
_, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
        method=BinningMethod(), maxrₑ=400., minrₑ=minrₑ
)

# 

using BSON: @save, @load

@load "fluxtest.bson" flux
plot(bins, flux)