include("utils/FittingUtils.jl")

bins = range(3, 10, 1000)

LPModel = XS_LampPostJohannsen(;
    K = .1,
    E = 6.4,
    a = 0.798,
    h = 19.,
    θ = 5.,
    α13 = 0.,
    ϵ3 = 5.,
) 

PLModel = PowerLaw(
    K = FitParam(50.),
    a = FitParam(2.)
)

flux = invokemodel!(bins, LPModel)
flux2 = invokemodel(bins, PLModel)
plot(bins, flux.+flux2)