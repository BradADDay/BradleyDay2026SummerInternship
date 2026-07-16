
using DataFrames
using CSV
using Dates
using ProgressBars
using Gradus

include("Deformations.jl")

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

function GetLineProfile(bins, a, h, θ, α13, ϵ3)

    # Position of the observer
    x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)

    # Instantiating the metric
    m = JohannsenMetric(M=1., a=a, α13=α13, ϵ3=ϵ3)

    # Disk
    d = ThinDisc(0., Inf)

    # Setting up the model and emissivity profile
    model = LampPostModel(h = h)
    profile = emissivity_profile(m, d, model)

    # Computing the line profile
    _, flux = lineprofile(m, x, d, profile; verbose=false, bins=bins/E, 
            method=TransferFunctionMethod()
    )

    return flux

end

function FixISCO(ϵ3, α13, a)

    while is_no_isco(ϵ3, α13, a)

        α13 += 0.1

    end

    return α13

end

function Generate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)
    # Looping through the parameter space, generating and saving spectra
    
    for a in as
        j=1
        df = DataFrame()
        for h in hs
            for θ in θs
                for α13 in α13s
                    strα13 = copy(α13)
                    for ϵ3 in ϵ3s

                        α13 = strα13

                        combination = "$h, $θ, $strα13, $ϵ3"
                        update(pbar)

                        try
                            
                            if (ϵ3 < Constraints(a))

                                ϵ3 = Constraints(a)

                                open(joinpath(OUTDIR, "log.txt"), "a") do io
                                    write(io, "$(now()): WARNING: Combination ($a, $h, $θ, $α13, $ϵ3) out of bounds!\n")
                                end

                            end
                            if (α13 < Constraints(a))

                                # Setting α13 to the lowest possible value
                                α13 = Constraints(a)

                                open(joinpath(OUTDIR, "log.txt"), "a") do io
                                    write(io, "$(now()): WARNING: Combination ($a, $h, $θ, $α13, $ϵ3) out of bounds!\n")
                                end

                            end
                            if is_no_isco(ϵ3, α13, a)

                                # Shifting α13 up until valid ISCO
                                α13 = FixISCO(ϵ3, α13, a)

                                open(joinpath(OUTDIR, "log.txt"), "a") do io
                                    write(io, "$(now()): WARNING: Combination ($a, $h, $θ, $α13, $ϵ3) has no ISCO! Setting α13=$α13\n")
                                end

                            end

                            open(joinpath(OUTDIR, "test.txt"), "a") do io
                                write(io, "$a, $h, $θ, $α13, $ϵ3\n")
                            end

                            # Calculating the spectrum and storing it in df
                            flux = GetLineProfile(bins, a, h, θ, α13, ϵ3)

                        catch err

                            # If the parameter combination fails, noting this in a log file
                            open(joinpath(OUTDIR, "log.txt"), "a") do io
                                write(io, "$(now()): Combination ($a, $h, $θ, $α13, $ϵ3) failed!\n")
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

function GetVars(path::String, as, hs, θs, α13s, ϵ3s)

    for i in [as, hs, θs, α13s, ϵ3s]
        open(joinpath(path, "vars.txt"), "a") do io
            write(io, "$(collect(i)),\n")
        end
    end

end

# Defining the parameter space
# as   = range(-0.998, 0.998, 12)
# hs   = range( 3  , 15.   , 9)
# θs   = range( 5.   , 85.   , 9)
# α13s = range(-8., 10., 10)
# ϵ3s  = range(-8., 10., 10)

as = [-0.998, 0., 0.998]
hs = [3., 9., 15.]
θs = [5., 45., 85.]
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

OUTDIR = "tabledata4/"

try
    mkdir(OUTDIR)
catch
end

GetVars(OUTDIR, as, hs, θs, α13s, ϵ3s)

# Using 1000 bins for high resolution to reduce the effects of interpolation
bins = collect(range(0., 3., 1000))

pbar = ProgressBar(total=length(hs)*length(as)*length(θs)*length(α13s)*length(ϵ3s))

MultiGenerate(as, pbar, bins, OUTDIR, hs, θs, α13s, ϵ3s)

# cp(OUTDIR, "/home/brad/OneDrive/$OUTDIR", force=true)