# =======================================================================
# Setup
# =======================================================================

# Metric Defaults
const M = 1.0     # Mass
const a = 0.998   # Spin

# Perturbation Defaults
const α13 = 0.
const α22 = 0.
const α52 = 0.
const ϵ3  = 0.

# BH Defaults
const θ = 60. # Inclination, degrees
const h = 10. # Corona height

# Dictionary for easy access and modification of the model parameters
defaultSetupDict = Dict((["θ", θ], ["α13", α13], ["M", M], ["α22", α22], ["ϵ3", ϵ3], ["a", a], ["h", h], ["α52", α52]))