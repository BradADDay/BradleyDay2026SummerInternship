begin
	"""
	Get a file of all combinations that fail
	"""

	using Gradus: JohannsenMetric, SVector, isco, ThinDisc, LampPostModel, emissivity_profile, lineprofile, TransferFunctionMethod, BinningMethod, transferfunctions
	# using OrdinaryDiffEqSDIRK: TRBDF2
	using ProgressBars

	include("Code/utils/Deformations.jl")

	function testing(as, α13s, ϵ3s, bins, method, maxrₑ, numrₑ, x, model, pbar)
		for a in as
			for α13 in α13s, ϵ3 in ϵ3s
				update(pbar)
				try

					if !QuickIsValid(ϵ3, α13, a)
						ϵ3, α13 = FindNearestSafePoint([ϵ3, α13], a)
					end

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
end

begin
	"""
	Generate JSON of fixed parameter combinations
	"""

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

	FixedParamsJSON(;fails=fails)
end

begin
	"""
	Testing Gradus
	"""

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

end

begin

	using Gradus: JohannsenMetric

	function draw_horizon(p, m)
		rs, θs = event_horizon(m, resolution = 200)
		radius = rs

		x = @. radius * sin(θs)
		y = @. radius * cos(θs)
		plot!(p, x, y, label = nothing)
	end

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

	p = plot(aspect_ratio = 1)
	for (a, α13, ϵ3) in fails
		m = JohannsenMetric(M = 1.0, a=a, ϵ3=ϵ3, α13=α13)
		draw_horizon(p, m)
	end
	p

end