from heasp import *
import numpy as np
import os
import pandas as pd

tbl = table()

def getColumn(a, h, theta, alpha, epsilon):
    return f"{a}, {h}, {theta}, {alpha}, {epsilon}"

# Primary header

tbl.setModelName("JhnsnLmPst")
tbl.setModelUnits(" ")
tbl.setisRedshift(False)
tbl.setisAdditive(True)
tbl.setisError(False)
tbl.setisEscale(True)

# Energies

tbl.setEnergies(np.linspace(0.2, 1.5, 1001))

# Parameters

tbl.setNumIntParams(5)
tbl.setNumAddParams(0)

spins = [-0.998, -0.6653333333333333, -0.33266666666666667, 0.0, 0.33266666666666667, 0.6653333333333333, 0.998]
spin = tableParameter("SPIN", 0, 0.998, 0.001, -0.998, -0.998, 0.998, 0.998)
spin.setTabulatedValues(spins)

tbl.pushParameter(spin)

heights = [1.5, 4.666666666666667, 7.833333333333333, 11.0, 14.166666666666666, 17.333333333333332, 20.5, 23.666666666666668, 26.833333333333332, 30.0]
height = tableParameter("HEIGHT", 0, 10.0, 0.1, 1.5, 1.5, 30.0, 30.0) 
height.setTabulatedValues(heights)

tbl.pushParameter(height)

inclinations = [5.0, 13.88888888888889, 22.77777777777778, 31.666666666666668, 40.55555555555556, 49.44444444444444, 58.333333333333336, 67.22222222222223, 76.11111111111111, 85.0]
inclination = tableParameter("INCLINATION", 0, 60.0, 0.1, 5.0, 5.0, 85.0, 85.0) 
inclination.setTabulatedValues(inclinations)

tbl.pushParameter(inclination)

alphas = [0.0, 8.333333333333334, 16.666666666666668, 25.0, 33.333333333333336, 41.666666666666664, 50.0]
alpha13 = tableParameter("ALPHA13", 0, 0.0, 0.1, 0.0, 0.0, 50.0, 50.0) 
alpha13.setTabulatedValues(alphas)

tbl.pushParameter(alpha13)

epsilons = [0.0, 3.3333333333333335, 6.666666666666667, 10.0, 13.333333333333334, 16.666666666666668, 20.0, 23.333333333333332, 26.666666666666668, 30.0]
epsilon3 = tableParameter("EPSILON3", 0, 0.0, 0.1, 0.0, 0.0, 30.0, 30.0)
epsilon3.setTabulatedValues(epsilons)

tbl.pushParameter(epsilon3)

df = pd.read_csv("data/spectra.csv")

i=0

for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:
                    try:
                        flux = df[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
                    except:
                        flux = np.zeros(1000)
                        i+=1
                    
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

tablefile = "test.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

print(i)