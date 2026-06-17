using Gradus

function computeLineProfile(m, x, height, bins = range(0.0, 1.5, 180))

    # Disk
    d = ThinDisc(10., Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = height)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=true, bins=bins, 
            method=TransferFunctionMethod(), minrₑ=10.)

    #= profile(r) = r^-2

    _, flux = lineprofile(bins, profile, m, x, d; verbose=true,
            method=TransferFunctionMethod(), minrₑ=10.) =#

    return flux, Gradus.isco(m), profile
end

m = KerrMetric(1., 0.0)
jm = JohannsenMetric(1., 0.998, 0., 0., 0., 0.)

x = SVector(0., 10000., deg2rad(60), 0.)

flux = computeLineProfile(m, x, 50.)[1]
jflux = computeLineProfile(jm, x, 50.)[1]

bins = range(0, 1.5, 180)

plot(bins, flux, label="Kerr")
plot!(bins, jflux, label="Johannsen")