"""
quick script to find the model parameters that failed
"""

f = open("data/log.txt", "r")

include("ParameterVariations.jl")

bins = collect(range(0.2, 1.5, 1000))

test = fill(zeros(5), 3500)

i=1

for line in readlines(f)
    a, h, θ, α13, ϵ3 = [parse(Float64, String(x)) for x in split(line[39:end-9], ", ")]
    test[i] = [a, h, θ, α13, ϵ3]

    if i % 50 == 0
        println("$i: $(test[i])")
        setup = Dict((
            ["θ", θ], 
            ["α13", -1.], 
            ["M", 1.], 
            ["α22", 0.], 
            ["ϵ3", -1.], 
            ["a", a], 
            ["h", 3.], 
            ["α52", 0.]
        ))

        # Calculating the spectrum and storing it in df
        flux = JohannsenParamVar(setup, bins, ComputeLineProfile)
    end
    i+=1
end

test = mapreduce(permutedims, vcat, test)


for i in 1:5
    println(unique(test[:,i]))
end

"""
a   : (all) [-0.6653333333333333, -0.33266666666666667, 3.266666666666667e15, 0.0, 0.33266666666666667, 2.66666666666667e14, 0.6653333333333333, 5.3333333333333e13]
h   :       [1.5]
θ   : (all) [5.0, 13.88888888888889, 22.77777777777778, 31.666666666666668, 40.55555555555556, 49.44444444444444, 58.333333333333336, 67.22222222222223, 76.11111111111111, 85.0]
α13 : (all) [0.0, 8.333333333333334, 16.666666666666668, 25.0, 33.333333333333336, 41.666666666666664, 50.0]
ϵ3  : (all) [0.0, 3.3333333333333335, 6.666666666666667, 10.0, 13.333333333333334, 16.666666666666668, 20.0, 23.333333333333332, 26.666666666666668, 30.0]
"""