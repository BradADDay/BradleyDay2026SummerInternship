using Dates
using ProgressBars
using JSON3
using Gradus: integrate_lineprofile

include("utils/Deformations.jl")

println("$(Threads.nthreads()) threads in use")

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

# Reading the JSON file containing the corrected coordinates in parameter space
const jsonFile = JSON3.read("Code/utils/FixedCoords.json")

"""
Generate CSVs of spectra to compile a table model for the Johannsen metric with a lamppost corona
"""
function Generate(as, α13s, ϵ3s, θs, hs, pbar, OUTDIR; bins=range(0, 3, 1000), json=jsonFile, maxrₑ=400.)

    # Dictionary for storing spectra to avoid unnecessary calculations of repeats
    spectra = Dict()

    for a in as
        
        # Variable to store the column index when appending to df
        j=1

        try
            df = DataFrame(CSV.File(joinpath(OUTDIR, "$(a).csv")))
        catch
            df = DataFrame()
        end

        FIXED = json["$(round(a; digits=3))"]

        for α13 in α13s, ϵ3 in ϵ3s, θ in θs

            # Storing the combination before the deformation parameters get overwritten
            combo = "$α13, $ϵ3, $θ"

            # Fixing the coordinates
            α13, ϵ3 = Float64.(FIXED["$α13, $ϵ3"])
            if α13 <= -5.7
                α13 = -5.7
            end

            # Setting up the metric, observer position, and disc
            m = JohannsenMetric(; M=1., a=a, α13=α13, ϵ3=ϵ3)
            x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)
            minrₑ = Gradus.isco(m)
            d = ThinDisc(minrₑ, Inf)

            # Pre-computing transfer functions
            tfs = transferfunctions(
                m, x, d;
                maxrₑ= maxrₑ,
                numrₑ=200, 
                minrₑ=minrₑ,
            )

            for h in hs

                update(pbar)

                # Storing the combination for error reporting
                combination = "$combo, $h"

                # Instantiating the flux
                flux = Array{Float64}(undef, length(bins))

                try

                    # Checking if there is a spectrum already computed
                    flux .= spectra[combination]

                catch
                    try

                        # Setting up the model and emissivity profile
                        model = LampPostModel(h = h)
                        profile = emissivity_profile(m, d, model)

                        # Computing the line profile
                        flux .= integrate_lineprofile(
                            profile, tfs, bins;
                            rmax=maxrₑ, 
                            rmin=minrₑ
                        )

                        # Storing the flux to avoid unnecessary computation
                        spectra[combination] = flux

                    catch err

                        # If the parameter combination fails, noting this in a log file
                        open(joinpath(OUTDIR, "$(a)_log.txt"), "a") do io
                            write(io, "$(now()): Combination ($combination) failed!\n")
                        end

                        # Noting the relevant error message
                        open(joinpath(OUTDIR, "$(a)_err.txt"), "a") do io
                            write(io, "$(now()): $err\n")
                        end

                        # Setting the flux to an array of zeros
                        flux = zeros(1000)

                        # Storing the flux to avoid unnecessary computations
                        spectra[combination] = flux
                    end
                end

                insertcols!(df, j, combination => flux)
                j+=1
                
            end

            # Saving to CSV periodically to avoid losing all data in the event of a crash etc.
            CSV.write(joinpath(OUTDIR, "$(a).csv"), df)

            df = DataFrame(CSV.File(joinpath(OUTDIR, "$(a).csv")))

            open(joinpath(OUTDIR, "$(round(a; digits=3))_spectra.json"), "w") do io
                JSON3.pretty(io, spectra)
            end
        
        end
    end
end

function MultiGenerate(as, α13s, ϵ3s, θs, hs, pbar, OUTDIR; bins=range(0, 3, 1000), json=jsonFile)

    """
    Dispatch `Generate(...)` to multiple threads
    """

    chunks = Iterators.partition(as, cld(length(as), Threads.nthreads()))
    tasks = map(chunks) do chunk
        Threads.@spawn Generate(chunk, α13s, ϵ3s, θs, hs, pbar, OUTDIR; bins=bins, json=json)
    end
    fetch.(tasks)

end

function GetVars(path::String, as, α13s, ϵ3s, θs, hs)

    """
    Generate a text file containing the parameters used, this is used to later construct the table model
    """

    for i in [as, hs, θs, α13s, ϵ3s]
        open(joinpath(path, "vars.txt"), "a") do io
            write(io, "$(collect(i)),\n")
        end
    end

end

# Setting up the parameter space to compute for
as   = range(0, 0.998, 10)
hs   = range( 3.  , 19.   , 9)
θs   = range( 5.   , 85.   , 9)
α13s = [10.]
ϵ3s  = reverse(range(-8., 10., 10))

# Output directory
OUTDIR = "tabledataFinal5/"

try
    mkdir(OUTDIR)
catch
end

# Generating the variables file
GetVars(OUTDIR, as, α13s, ϵ3s, θs, hs)

# Setting up a progress bar
pbar = ProgressBar(total=length(hs)*length(θs)*length(ϵ3s)*length(α13s)*length(as))

# Generating the CSVs
MultiGenerate(as, α13s, ϵ3s, θs, hs, pbar, OUTDIR)

cp(OUTDIR, "/home/brad/OneDrive/tabledataFinal5", force=true)

