using SpectralFitting, Plots
using Random
Random.seed!(0)

x = collect(range(-5, 5, 200))

A_true = 3.0
pos_true = 1.3
sigma_true = 0.8
err_true = 0.2

y = @. A_true * exp(-(x - pos_true)^2 / (2 * sigma_true^2))

y_noisy = y .+ (0.2 * randn(length(y)))

data = InjectiveData(x, y_noisy; name = "example")

plot(data, markersize = 3)

model = GaussianLine(μ = FitParam(0.0))

plot(data.domain[1:end-1], invokemodel(data.domain, model))

prob = FittingProblem(model => data)

result = fit(prob, LevenbergMarquadt())

plot(data, markersize = 3)
plot!(result)

amps = range(50, 200, 50)
devs = range(0.5, 1.2, 50)

stats = [
    measure(ChiSquared(), result, [a, result.u[2], d])
    for d in devs, a in amps
]

# 1, 2, and 3 sigma contours
stdev = std(stats)
levels = [1stdev, 2stdev, 3stdev]
contour(
    amps,
    devs,
    stats .- sum(result.stats),
    levels = levels,
    xlabel = "K",
    ylabel = "σ"
)
scatter!([result.u[1]], [result.u[3]])