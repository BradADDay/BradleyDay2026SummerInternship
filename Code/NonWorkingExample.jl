using Gradus

# Setup
M = 1.
a = 0.998
h = 10.
θ = 60.
α13 = -0.4
ϵ3 = -0.4

bins = range(3, 10, 1000)

# Position of the observer
x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

# Instantiating the metric
m = JohannsenMetric(M=M, a=a, α13=α13, ϵ3=ϵ3)

# Disk
d = ThinDisc(0., Inf)

# Setting up the model and emissivity profile
model = LampPostModel(h=h)
profile = emissivity_profile(m, d, model)

# Computing the line profile
_, flux = lineprofile(m, x, d, profile; verbose=false, bins=bins, method=TransferFunctionMethod())
