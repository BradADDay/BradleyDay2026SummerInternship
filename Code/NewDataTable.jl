using Dates
using ProgressBars
using JSON3

include("utils/Deformations.jl")

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

global const jsonFile = JSON3.read("Code/utils/FixedCoords.json")

function GetLineProfile(bins, tfs, h, minrₑ)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = h)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(
        profile, tfs; verbose=false, bins=bins, 
        method=TransferFunctionMethod(), maxrₑ=400, 
        numrₑ=50, minrₑ=minrₑ
    )

    return flux

end

function Generate(as, α13s, ϵ3s, θs, hs, pbar, bins, OUTDIR, json)

    for a in as
        
        j=1
        df = DataFrame()

        for α13 in α13s, ϵ3 in ϵ3s, θs in θs

            println("$α13, $ϵ3")

            if !QuickIsValid(ϵ3, α13, a)
                ϵ3, α13 = convert.(Float64, json["$ϵ3, $α13"])
            end

            m = JohannsenMetric(; M=1., a=a, α13=α13, ϵ3=ϵ3)
            x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)
            d = ThinDisc(Gradus.isco(m), Inf)

            # tfs = transferfunctions(m, x, d)
        
        end
    end

end

as   = range(-0.998, 0.998, 13)
hs   = range( 3  , 15.   , 8)
θs   = range( 5.   , 85.   , 8)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

bins = collect(range(0., 3., 1000))
OUTDIR = "FinalTableData/"
pbar = ProgressBar(total=length(hs)*length(as)*length(θs)*length(α13s)*length(ϵ3s))

Generate(as, hs, θs, α13s, ϵ3s, pbar, bins, OUTDIR, jsonFile)

# for h in hs
#     flux = GetLineProfile(bins, tfs, h, minrₑ)
# end