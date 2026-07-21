using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions
using Plots
using OrdinaryDiffEqBDF: FBDF

# Defining parameters
bins = range(0, 3, 1000)

a = -0.998
h = 3.0
θ = 45.0
# α13 = -0.99
# ϵ3 = 1.0
α13 = 0.
ϵ3 = 0.

method = TransferFunctionMethod()
maxrₑ = 400.
numrₑ = 200
maxiters = 1e7

println("position")
x = SVector(
        0.0, 
        10000.0, 
        deg2rad(θ), 
        0.0
)

println("metric")
m = JohannsenMetric(
        M=1., 
        a=a, 
        α13=α13, 
        ϵ3=ϵ3
)

println("disk")
minrₑ = isco(m)
d = ThinDisc(minrₑ, Inf)

tfs = transferfunctions(m, x, d; minrₑ=minrₑ, maxrₑ=maxrₑ, numrₑ=numrₑ)

println("lamppost model")
model = LampPostModel(h=h)
profile = emissivity_profile(
        m, d, model; 
        verbose=true,
        maxiters=maxiters
)

println("line profile")
_, flux = lineprofile(
        m, x, d; 
        verbose=true, 
        bins=bins, 
        method=method, 
        maxrₑ=maxrₑ, 
        numrₑ=numrₑ, 
        minrₑ=minrₑ,
        maxiters=maxiters
)

plot(bins, flux)