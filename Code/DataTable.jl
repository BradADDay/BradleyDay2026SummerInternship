
<<<<<<< HEAD
using DataFrames
using CSV
using Dates

"""
Generate a csv of line profiles by varying the 5 parameters for the Johannsen metric:
    a, h, θ, α13, ϵ3
"""

=======
using FITSFiles
using Random
using DataFrames
using CSV
using BenchmarkTools
using Dates

>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
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

<<<<<<< HEAD
# Defining the parameter space
hs   = range( 1.5  , 30.   , 10)
as   = range(-0.998,  0.998, 7)
θs   = range( 5.   , 85.   , 10)
α13s = range(0.    , 50.   , 7)
ϵ3s  = range(0.    , 30.   , 10)

# Output file
file = "output/spectra.csv"
df = DataFrame()
i=1

# Using 1000 bins for high resolution to reduce the effects of interpolation
bins = collect(range(0.2, 1.5, 1000))

# Looping through the parameter space, generating and saving spectra
=======
num = 2

hs   = range( 1.5  , 30.   , 10)
as   = range(-0.998,  0.998, 7)
θs   = range( 5.   , 85.   , 10)
α13s = range(0.   , 50.   , 7)
ϵ3s  = range(0.   , 30.   , 10)

hs   = range( 1.5  , 30.   , 2)
as   = range(-0.998,  0.998, 2)
θs   = range( 5.   , 85.   , 2)
α13s = range(0.   , 50.   , 2)
ϵ3s  = range(0.   , 30.   , 2)

file1 = "output/spectra.csv"

df = DataFrame()

i=1

bins = collect(range(0.2, 1.5, 1000))

>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
for a in as
    for h in hs
        for θ in θs
            for α13 in α13s
                for ϵ3 in ϵ3s
                    combination = "$a, $h, $θ, $α13, $ϵ3"
                    try
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

<<<<<<< HEAD
                        # Calculating the spectrum and storing it in df
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
                        flux = JohannsenParamVar(setup, bins, ComputeLineProfile)
                        insertcols!(df, i, combination => flux)
                        i+=1
                    catch
<<<<<<< HEAD
                        # If the parameter combination fails, noting this in a log file
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
                        open("output/log.txt", "a") do io
                            write(io, "$(now()): Combination ($combination) failed!\n")
                        end
                    end
                end
            end
<<<<<<< HEAD
            # Saving to CSV periodically to avoid losing all data in the event of a crash etc.
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
            CSV.write(file1, df)
            df = DataFrame(CSV.File(file1))
        end
    end
end