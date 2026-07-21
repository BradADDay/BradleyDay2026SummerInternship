using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions
# using OrdinaryDiffEqSDIRK: TRBDF2
using ProgressBars

include("Code/utils/Deformations.jl")

function testing(as, α13s, ϵ3s, bins, method, maxrₑ, numrₑ, x, model, pbar)
	for a in as
		for α13 in α13s, ϵ3 in ϵ3s
			update(pbar)
			try
				# println()
				# Defining parameters
				
				# print("$a, $α13, $ϵ3 => ")

				if !QuickIsValid(ϵ3, α13, a)
						ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
				end

				# println("$a, $α13, $ϵ3")

				# println("Valid: $(IsValid(ϵ3, α13, a))")

				m = JohannsenMetric(
					M=1., 
					a=a, 
					α13=α13, 
					ϵ3=ϵ3
				)

				minrₑ = isco(m)
				d = ThinDisc(minrₑ, Inf)

				profile = emissivity_profile(
					m, d, model
				)

				_, flux = lineprofile(
					m, x, d, profile; 
					bins=bins, 
					method=method, 
					maxrₑ=maxrₑ, 
					numrₑ=numrₑ, 
					minrₑ=minrₑ
				)

			catch err
				print("fail")
				open("outBin/testing.txt", "a") do io
					write(io, "$a, $α13, $ϵ3: $err\n")
				end
			end
		end
	end
end

as   = range(0, 0.998, 10)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

bins = range(0, 3, 1000)
h = 3.0
θ = 45.0

method = TransferFunctionMethod()
maxrₑ = 400.
numrₑ = 100

x = SVector(
	0.0, 
	10000.0, 
	deg2rad(θ), 
	0.0
)

model = LampPostModel(h=h)

pbar = ProgressBar(total = 1000)

chunks = Iterators.partition(as, cld(length(as), Threads.nthreads()))
tasks = map(chunks) do chunk
	Threads.@spawn testing(as, α13s, ϵ3s, bins, method, maxrₑ, numrₑ, x, model, pbar)
end
fetch.(tasks)