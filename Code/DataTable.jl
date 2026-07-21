using Dates
using ProgressBars
using JSON3

include("utils/Deformations.jl")

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

global const jsonFile = JSON3.read("Code/utils/FixedCoords.json")

function GetLineProfile(bins, a, h, θ, α13, ϵ3)

    x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

    # Instantiating the metric
    m = JohannsenMetric(M=1., a=a, α13=α13, ϵ3=ϵ3)

    minrₑ = Gradus.isco(m)

    d = ThinDisc(minrₑ, Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = h)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=false, bins=bins, 
            method=TransferFunctionMethod(), maxrₑ=400, numrₑ=50, minrₑ=minrₑ
    )

    return flux

end

function Generate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)
    # Looping through the parameter space, generating and saving spectra

    for a in as
    
        FixedCoords = jsonFile["$(round(a; digits=3))"]
        j=1
        df = DataFrame()

        for h in hs
            for θ in θs
                for α13 in α13s
                    strα13 = copy(α13)
                    for ϵ3 in ϵ3s

                        α13 = strα13
                        combination = "$h, $θ, $strα13, $ϵ3"

                        try

                            if !QuickIsValid(ϵ3, α13, a)
                                ϵ3, α13 = convert.(Float64, FixedCoords["$ϵ3, $α13"])
                            end

                            # Calculating the spectrum and storing it in df
                            flux = GetLineProfile(bins, a, h, θ, α13, ϵ3)
                            
                        catch err
                            # If the parameter combination fails, noting this in a log file
                            open(joinpath(OUTDIR, "log.txt"), "a") do io
                                write(io, "$(now()): Combination ($combination) failed!\n")
                            end

                            open(joinpath(OUTDIR, "err.txt"), "a") do io
                                write(io, "$(now()): $err\n")
                            end

                            flux = zeros(1000)
                        end

                        insertcols!(df, j, combination => flux)
                        j+=1
                        update(pbar)

                    end
                end
            end
            # Saving to CSV periodically to avoid losing all data in the event of a crash etc.
            CSV.write(joinpath(OUTDIR, "$(a).csv"), df)
        end
    end
end

function MultiGenerate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)

    chunks = Iterators.partition(as, cld(length(as), Threads.nthreads()))
    tasks = map(chunks) do chunk
        Threads.@spawn Generate(chunk, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)
    end
    fetch.(tasks)

end

function GetVars(path::String, as, hs, θs, α13s, ϵ3s)

    for i in [as, hs, θs, α13s, ϵ3s]
        open(joinpath(path, "vars.txt"), "a") do io
            write(io, "$(collect(i)),\n")
        end
    end

end

# Defining the parameter space
as   = range(-0.998, 0.998, 13)
hs   = range( 3  , 15.   , 8)
θs   = range( 5.   , 85.   , 8)
as   = range(-0.998, 0.998, 13)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

OUTDIR = "FinalTableData/"

try
    mkdir(OUTDIR)
catch
end

GetVars(OUTDIR, as, hs, θs, α13s, ϵ3s)

# Using 1000 bins for high resolution to reduce the effects of interpolation
bins = collect(range(0., 3., 1000))

pbar = ProgressBar(total=length(hs)*length(as)*length(θs)*length(α13s)*length(ϵ3s))

MultiGenerate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)

cp(OUTDIR, "/home/brad/OneDrive/$OUTDIR", force=true)