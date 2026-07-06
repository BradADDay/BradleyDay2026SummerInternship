
using DataFrames
using CSV
using Dates
using ProgressBars

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

include("ParameterVariations.jl")

"""
h   ∈ [1.5, 30]
θ   ∈ [5, 85]
a   ∈ [-0.998, 0.998]
α13 ∈ [-(1+sqrt(1 - a^2))^3, 50]
ϵ3  ∈ [-(1+sqrt(1 - a^2))^3, 30]

approx 1.4179 per parameter combo
39 hours for 10 iterations of each parameter
"""

# Defining the parameter space
as   = range(0, 0.998, 5)
hs   = range( 3  , 15.   , 7)
θs   = range( 5.   , 85.   , 7)
α13s = [-0.4, 0., 2., 4., 6., 8., 10.]
ϵ3s  = [-0.4, 0., 2., 4., 6., 8., 10.]

for i in [as, hs, θs, α13s, ϵ3s]
    open("vars.txt", "a") do io
        write(io, "$(collect(i)),\n")
    end
end

OUTDIR = "tabledata/"

# mkdir(OUTDIR)

# Output file
df = DataFrame()
i=1
j=1

# Using 1000 bins for high resolution to reduce the effects of interpolation
bins = collect(range(0., 3., 1000))

pbar = ProgressBar(total=length(hs)*length(as)*length(θs)*length(α13s)*length(ϵ3s))

# Looping through the parameter space, generating and saving spectra
for a in as
    for h in hs
        for θ in θs
            for α13 in α13s
                for ϵ3 in ϵ3s
                    combination = "$a, $h, $θ, $α13, $ϵ3"
                    try

                        update(pbar)

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
                        insertcols!(df, i, combination => flux)
                        i+=1
                    catch err
                        # If the parameter combination fails, noting this in a log file
                        open(joinpath(OUTDIR, "log.txt"), "a") do io
                            write(io, "$(now()): Combination ($combination) failed!\n")
                        end
                        open(joinpath(OUTDIR, "err.txt"), "a") do io
                            write(io, "$(now()): $err\n")
                        end
                    end
                end
            end
        end
    end
    # Saving to CSV periodically to avoid losing all data in the event of a crash etc.
    CSV.write(joinpath(OUTDIR, "$j"), df)
    j+=1
    df = DataFrame()
    i=1
end