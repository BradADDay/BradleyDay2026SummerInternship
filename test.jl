using Plots

include("LampPostModelFit.jl")
include("ParameterVariations.jl")
include("Defaults.jl")

n = 0

setupDict = copy(defaultSetupDict)

fluxes = Array{Any}(nothing, 32)
configs = Array{Any}(nothing, 32)

i=1

for θ in [0.1, 89.9]
    for a in [0.05, 0.998]
        for h in [5, 30]
            for α13 in [0, 50]
                for ϵ3 in [0, 30]
                    setupDict["θ"] = θ
                    setupDict["a"] = a
                    setupDict["h"] = h
                    setupDict["α13"] = α13
                    setupDict["ϵ3"] = ϵ3
                    fluxes[i] = JohannsenParamVar(setupDict, range(0, 1.5, 180), ComputeLineProfile)
                    configs[i] = generateConfig(["θ", "a", "h", "α13", "ϵ3"], setupDict)
                    display(plot(range(0, 1.5, 180), fluxes[i], title=configs[i]))
                    i+=1
                end
            end
        end
    end
end