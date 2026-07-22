
# using OrdinaryDiffEqSDIRK: TRBDF2
using ProgressBars

include("Code/utils/Deformations.jl")
include("Code/utils/DeformUtils.jl")

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

for (i, (a, α13, ϵ3)) in enumerate(fails)
	
	hmp = PlotError(ϵ3, α13, a)
	savefig(hmp, "$i.png")

	close("all")
	GC.gc()
	
end

as   = range(0, 0.998, 10)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

##

using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions

fails = [
	[0.33266666666666667, -6.84, -5.96], 
	[0.6653333333333333, -4.8, -3.92], 
	[0.5544444444444444, -5.55, -3.94], 
	[0.6653333333333333, -4.71, -1.46], 
	[0.5544444444444444, -5.55, -3.93], 
	[0.7762222222222223, -4.08, -1.58], 
	[0.8871111111111111, -2.61, 0.03], 
	[0.998, -1.02, 6.09], 
	[0.7762222222222223, -4.23, 0.1], 
	[0.7762222222222223, -4.2, 2.1], 
	[0.7762222222222223, -3.92, 4.05], 
	[0.7762222222222223, -3.81, -1.72], 
	[0.7762222222222223, -4.0, 0.0], 
	[0.7762222222222223, -4.0, 2.0], 
	[0.7762222222222223, -3.9, 4.1]
]

include("Code/utils/DeformUtils.jl")

for (i, (a, α13, ϵ3)) in enumerate(fails)

	α13 += 0.7

	hmp = PlotError(ϵ3, α13, a)
	display(hmp)

	close("all")
	GC.gc()
	
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

	m = JohannsenMetric(
		M=1., 
		a=a, 
		α13=α13, 
		ϵ3=ϵ3
	)

	minrₑ = isco(m)
	d = ThinDisc(minrₑ, Inf)

	model = LampPostModel(h=h)
	profile = emissivity_profile(
		m, d, model
	)

	_, flux = lineprofile(
		m, x, d, profile; 
		bins=bins, 
		method=method, 
		maxrₑ=maxrₑ, 
		numrₑ=numrₑ, 
		minrₑ=minrₑ,
		zero_atol=1e-5
	)
end

##

as   = range(0, 0.998, 10)
α13s = range(-8., 10., 10)
ϵ3s  = range(-8., 10., 10)

fails = [
	[0.33266666666666667, -6.84, -5.96], 
	[0.6653333333333333, -4.8, -3.92], 
	[0.5544444444444444, -5.55, -3.94], 
	[0.6653333333333333, -4.71, -1.46], 
	[0.5544444444444444, -5.55, -3.93], 
	[0.7762222222222223, -4.08, -1.58], 
	[0.8871111111111111, -2.61, 0.03], 
	[0.998, -1.02, 6.09], 
	[0.7762222222222223, -4.23, 0.1], 
	[0.7762222222222223, -4.2, 2.1], 
	[0.7762222222222223, -3.92, 4.05], 
	[0.7762222222222223, -3.81, -1.72], 
	[0.7762222222222223, -4.0, 0.0], 
	[0.7762222222222223, -4.0, 2.0], 
	[0.7762222222222223, -3.9, 4.1]
]

include("Code/utils/DeformUtils.jl")

dict = Dict()
f=0

for a in as
	dict["$(round(a; digits=3))"] = Dict()
	for α13 in α13s, ϵ3 in ϵ3s

		iα13 = α13
		iϵ3 = ϵ3

		if !QuickIsValid(ϵ3, α13, a)
			ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
		end

		if [a, α13, ϵ3] in fails
			α13 += 0.6
			f += 1
		end

		dict["$(round(a; digits=3))"]["$iα13, $iϵ3"] = (α13, ϵ3)

		println("$a, $iα13, $iϵ3 => $a, $α13, $ϵ3")

	end
end

dict["numFails"] = f

println(f)

using JSON3

open("Code/utils/FixedCoords.json", "w") do f
	JSON3.pretty(f, dict)
	println(f)
end

##

using Gradus
using Plots

a=0.998
α13=0.
ϵ3=0.
h=10.
θ=60.
maxrₑ = 400.

bins = range(0, 3, 1000)

# Pre computing transfer functions
m = JohannsenMetric(; M=1., a=a, α13=α13, ϵ3=ϵ3)
x = SVector(0.0, 10000.0, deg2rad(θ), 0.0)
minrₑ = Gradus.isco(m)
d = ThinDisc(minrₑ, Inf)

tfs = transferfunctions(
	m, x, d;
	maxrₑ= maxrₑ,
	numrₑ=200, 
	minrₑ=minrₑ,
)

model = LampPostModel(h = h)
profile = emissivity_profile(m, d, model)

# Computing the line profile
flux = integrate_lineprofile(
	profile, tfs, bins;
	rmax=maxrₑ, 
	rmin=minrₑ
)

plot(bins, flux)
xlims!(0.1, 1.5)

##

using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions

fails = [
	[0.33266666666666667, -6.84, -5.96], 
	[0.6653333333333333, -4.8, -3.92], 
	[0.5544444444444444, -5.55, -3.94], 
	[0.6653333333333333, -4.71, -1.46], 
	[0.5544444444444444, -5.55, -3.93], 
	[0.7762222222222223, -4.08, -1.58], 
	[0.8871111111111111, -2.61, 0.03], 
	[0.998, -1.02, 6.09], 
	[0.7762222222222223, -4.23, 0.1], 
	[0.7762222222222223, -4.2, 2.1], 
	[0.7762222222222223, -3.92, 4.05], 
	[0.7762222222222223, -3.81, -1.72], 
	[0.7762222222222223, -4.0, 0.0], 
	[0.7762222222222223, -4.0, 2.0], 
	[0.7762222222222223, -3.9, 4.1]
]

a, α13, ϵ3 = fails[1]

bins = range(0, 3, 1000)
h = 3.0
θ = 45.0

method = TransferFunctionMethod()
maxrₑ = 400.
numrₑ = 200

x = SVector(
	0.0, 
	10000.0, 
	deg2rad(θ), 
	0.0
)

m = JohannsenMetric(
	M=1., 
	a=a, 
	α13=α13, 
	ϵ3=ϵ3
)

minrₑ = isco(m)
d = ThinDisc(minrₑ, Inf)

model = LampPostModel(h=h)
profile = emissivity_profile(
	m, d, model
)

_, flux = lineprofile(
	m, x, d, profile; 
	bins=bins, 
	method=method, 
	maxrₑ=maxrₑ, 
	numrₑ=numrₑ, 
	minrₑ=minrₑ,
	verbose=true,
	zero_atol=1e-2
)