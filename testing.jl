using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions
using Plots

include("Code/utils/Deformations.jl")

for i in 1:10
        # Defining parameters
        bins = range(0, 3, 1000)

        a = rand(0:0.001:0.998)
        h = 3.0
        θ = 45.0
        α13 = rand(-8:0.1:10)
        ϵ3 = rand(-8:0.1:10)

        println("$a, $α13, $ϵ3")

        if !QuickIsValid(ϵ3, α13, a)
                ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
        end

        method = TransferFunctionMethod()
        maxrₑ = 400.
        numrₑ = 200

        println("position")
        x = SVector(
        0.0, 
        1000.0, 
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

        println("lamppost model")
        model = LampPostModel(h=h)
        profile = emissivity_profile(
        m, d, model; 
        verbose=true
        )

        println("line profile")
        _, flux = lineprofile(
        m, x, d; 
        verbose=true, 
        bins=bins, 
        method=method, 
        maxrₑ=maxrₑ, 
        numrₑ=numrₑ, 
        minrₑ=minrₑ
        )

        plot(bins, flux)
end