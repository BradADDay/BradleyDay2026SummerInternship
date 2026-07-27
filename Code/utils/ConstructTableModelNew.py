from heasp import *
import numpy as np
import os
import pandas as pd

"""
Convert the csv generated from DataTable.jl into a FITS file following the XSPEC
table model specifications. 

See https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/ogip_92_009/ogip_92_009.pdf
for more information.

This makes use of the python wrapper for the heasp C module.

See https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/headas/heasp/node1.html
for more information.
"""

0.8165454545454546
0.8165454545454546

spins, heights, inclinations, alphas, epsilons = [
    [0.0, 0.11088888888888888, 0.22177777777777777, 0.33266666666666667, 0.44355555555555554, 0.5544444444444444, 0.6653333333333333, 0.7762222222222223, 0.8871111111111111, 0.998],
    [1.0, 2.5555555555555554, 4.111111111111111, 5.666666666666667, 7.222222222222222, 8.777777777777779, 10.333333333333334, 11.88888888888889, 13.444444444444445, 15.0],
    [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0],
    [10.0, 8.0, 6.0, 4.0, 2.0, 0.0, -2.0, -4.0, -6.0, -8.0],
    [10.0, 8.0, 6.0, 4.0, 2.0, 0.0, -2.0, -4.0, -6.0, -8.0]
]

tbl = table()

# Utility function to generate a string formatted like the columns
def getColumn(h, theta, alpha, epsilon):
    return f"{alpha}, {epsilon}, {h}, {theta}"

# ===============================================================
# Primary header unit
# ===============================================================

tbl.setModelName("JhnsnLmPst")
tbl.setModelUnits(" ")
tbl.setisRedshift(False)
tbl.setisAdditive(True)
tbl.setisError(False)
tbl.setisEscale(True)

# ===============================================================
# Defining the energies
# ===============================================================

tbl.setEnergies(np.linspace(0., 3., 1001))

# ===============================================================
# Defining the parameters
# ===============================================================

tbl.setNumIntParams(5)
tbl.setNumAddParams(0)

# ---------------------------------------------------------------
# Spin
# ---------------------------------------------------------------

spin = tableParameter("SPIN", 0, 0.998, 0.001, min(spins), min(spins), max(spins), max(spins))
spin.setTabulatedValues(spins)

tbl.pushParameter(spin)

# ---------------------------------------------------------------
# Height
# ---------------------------------------------------------------

height = tableParameter("HEIGHT", 0, 10.0, 0.1, min(heights), min(heights), max(heights), max(heights)) 
height.setTabulatedValues(heights)

tbl.pushParameter(height)

# ---------------------------------------------------------------
# Inclination
# ---------------------------------------------------------------

inclination = tableParameter("INCLINATION", 0, 60.0, 0.1, min(inclinations), min(inclinations), max(inclinations), max(inclinations)) 
inclination.setTabulatedValues(inclinations)

tbl.pushParameter(inclination)

# ---------------------------------------------------------------
# Alpha13
# ---------------------------------------------------------------

alpha13 = tableParameter("ALPHA13", 0, 0.0, 0.1, min(alphas), min(alphas), max(alphas), max(alphas)) 
alpha13.setTabulatedValues(alphas)

tbl.pushParameter(alpha13)

# ---------------------------------------------------------------
# Epsilon3
# ---------------------------------------------------------------

epsilon3 = tableParameter("EPSILON3", 0, 0.0, 0.1, min(epsilons), min(epsilons), max(epsilons), max(epsilons))
epsilon3.setTabulatedValues(epsilons)

tbl.pushParameter(epsilon3)

# ===============================================================
# Reading the data to store in the table model
# ===============================================================

# Reading the CSV
print("Loading CSVs...")
DIR = "tabledataFinal4/"

files = {}

for i in os.listdir(DIR):
    if i[-4:] == ".csv":
        print(i)
        files[i[:-4]] = pd.read_csv(DIR + i)

i=0
j=0
k=0

# Looping through the parameters in the same order as for generation
print("Parameter loop...")
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:

                    j+=1

                    try:
                        flux = files[str(a)][getColumn(h, theta, alpha, epsilon)].to_numpy()
                    except:
                        flux = None
                    
                    if flux is None:
                        flux = np.zeros(1000)
                        i+=1

                    # Storing the spectrum and pushing it to the table alongside its parameter combination
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

# Saving the file
tablefile = DIR + "model.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

# Printing the number of failed spectra to compare against the logs
print(f"{100*(i/j)}% Failure")
print(f"{100*(k/j)}% NaN")
