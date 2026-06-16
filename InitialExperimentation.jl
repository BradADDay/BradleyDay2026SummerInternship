using Gradus
using Plots

# Function defining the cross section of the disk
function discCrossSection(x)
    centre = 8
    radius = 3

    if (x < centre - radius) || (radius + centre < x)
        zero(x)
    else 
        r = x - centre
        sqrt(radius^2 - r^2) + (0.5sin(3x))
    end
end

# Metric Parameters
M = 1.0
a = 0.7
α13 = 1.0
α22 = 0.0
α52 = 0.0
ϵ3 = 0.0
inclination = deg2rad(70)

# Instantiating the metric
m = JohannsenMetric(M, a, α13, α22, α52, ϵ3)

# Position of the observer
# (x[1], Distance, Inclination, x[4])
x = SVector(0.0, 1000.0, inclination, 0.0)

# Maximum affine time, parameterises the solution
# λ_max ~ 2*x[2]
λ_max = 2000.0

# Redshift point function
redshift = ConstPointFunctions.redshift(m, x)
redshiftGeometry = redshift ∘ ConstPointFunctions.filter_intersected()

# Disc
d = ThinDisc(0.0, 10.0)

# this function returns the impact parameter axes
α, β, image = rendergeodesics(
    m, x, d, λ_max, pf = redshiftGeometry,
    # image parameters
    image_width = 800, image_height = 600,
    # the "zoom" -- use the impact parameter axes
    αlims = (-20, 20), βlims = (-15, 15), verbose = true
)

heatmap(α, β, image, aspect_ratio = 1, c = :seismic)
