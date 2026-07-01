# =======================================================================
# Setup
# =======================================================================

# Metric Defaults
M = 1.0     # Mass
a = 0.998   # Spin

# Perturbation Defaults
α13 = 0.
α22 = 0.
α52 = 0.
ϵ3  = 0.

# BH Defaults
θ = 60. # Inclination, degrees
h = 10. # Corona height

# Dictionary for easy access and modification of the model parameters
defaultSetupDict = Dict((["θ", θ], ["α13", α13], ["M", M], ["α22", α22], ["ϵ3", ϵ3], ["a", a], ["h", h], ["α52", α52]))