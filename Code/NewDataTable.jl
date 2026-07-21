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

    

    return flux

end

function Generate(as, α13s, ϵ3s, θs, hs, pbar, bins, OUTDIR, json)

    # Dictionary for storing spectra to avoid unnecessary calculations of repeats
    spectra = Dict()

    for a in as
        
        # Variable to store the column index when appending to df
        j=1
        df = DataFrame()

        for α13 in α13s, ϵ3 in ϵ3s, θs in θs

            update(pbar)

            # Fixing the coordinates
            # TODO: change this to read from a table file
            if !QuickIsValid(ϵ3, α13, a)
                ϵ3, α13 = convert.(Float64, json["$ϵ3, $α13"])
            end

            # Pre computing transfer functions
            m = JohannsenMetric(; M=1., a=a, α13=α13, ϵ3=ϵ3)
            x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)
            minrₑ = Gradus.isco(m)
            d = ThinDisc(minrₑ, Inf)

            tfs = transferfunctions(m, x, d)

            for h in hs

                # Storing the combination for error reporting
                combination = "$a, $α13, $ϵ3, $θ, $h"

                try

                    # Checking if there is a spectrum already computed
                    flux = spectra[combination]

                catch
                    try

                        # Setting up the model and emissivity profile
                        model = LampPostModel(h = h)
                        profile = emissivity_profile(m, d, model)

                        # Computing the line profile
                        _, flux = lineprofile(
                            profile, tfs; verbose=false, bins=bins, 
                            method=TransferFunctionMethod(), maxrₑ=400, 
                            numrₑ=100, minrₑ=minrₑ
                        )

                        spectra[combination] = flux

                    catch err
                        # If the parameter combination fails, noting this in a log file
                        open(joinpath(OUTDIR, "log.txt"), "a") do io
                            write(io, "$(now()): Combination ($combination) failed!\n")
                        end

                        open(joinpath(OUTDIR, "err.txt"), "a") do io
                            write(io, "$(now()): $err\n")
                        end

                        flux = zeros(1000)

                        spectra[combination] = flux
                    end
                end

                insertcols!(df, j, combination => flux)
                j+=1
                
            end

            # Saving to CSV periodically to avoid losing all data in the event of a crash etc.
            CSV.write(joinpath(OUTDIR, "$(a).csv"), df)         
        
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