
using DataFrames
using CSV
using Dates
using ProgressBars

include("Deformations.jl")
include("ParameterVariations.jl")

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

function Generate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)
    # Looping through the parameter space, generating and saving spectra
    
    for a in as
        j=1
        df = DataFrame()
        for h in hs
            for θ in θs
                for α13 in α13s
                    for ϵ3 in ϵ3s
                        
                        combination = "$a, $h, $θ, $α13, $ϵ3"
                        update(pbar)
                        flux= fill(NaN, length(bins))

                        try
                            
                            if (ϵ3 < Constraints(a)) | (α13 < Constraints(a))

                                open(joinpath(OUTDIR, "log.txt"), "a") do io
                                    write(io, "$(now()): WARNING: Combination ($combination) out of bounds!\n")
                                end

                            elseif is_no_isco(ϵ3, α13, a)

                                open(joinpath(OUTDIR, "log.txt"), "a") do io
                                    write(io, "$(now()): WARNING: Combination ($combination) has no ISCO!\n")
                                end

                            else

                                setup = Dict((
                                    ["θ", θ], 
                                    ["α13", α13], 
                                    ["M", 1.], 
                                    ["α22", 0.], 
                                    ["ϵ3", ϵ3], 
                                    ["a", a], 
                                    ["h", h], 
                                    ["α52", 0.]
                                ))

                                # Calculating the spectrum and storing it in df
                                flux = JohannsenParamVar(setup, bins, ComputeLineProfile)

                            end

                        catch err

                            # If the parameter combination fails, noting this in a log file
                            open(joinpath(OUTDIR, "log.txt"), "a") do io
                                write(io, "$(now()): Combination ($combination) failed!\n")
                            end

                            open(joinpath(OUTDIR, "err.txt"), "a") do io
                                write(io, "$(now()): $err\n")
                            end

                        end

                        insertcols!(df, j, combination => flux)
                        j+=1
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

# Defining the parameter space
as   = range(-0.998, 0.998, 12)
hs   = range( 3  , 15.   , 8)
θs   = range( 5.   , 85.   , 8)
α13s = range(-8., 10., 12)
ϵ3s  = range(-8., 10., 12)

for i in [as, hs, θs, α13s, ϵ3s]
    open("vars.txt", "a") do io
        write(io, "$(collect(i)),\n")
    end
end

OUTDIR = "tabledata/"

# Using 1000 bins for high resolution to reduce the effects of interpolation
bins = collect(range(0., 3., 1000))

pbar = ProgressBar(total=length(hs)*length(as)*length(θs)*length(α13s)*length(ϵ3s))

MultiGenerate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)