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

# The values each parameter was varied between
spins, heights, inclinations, alphas, epsilons = [
    [0.0, 0.11088888888888888, 0.22177777777777777, 0.33266666666666667, 0.44355555555555554, 0.5544444444444444, 0.6653333333333333, 0.7762222222222223, 0.8871111111111111, 0.998],
    [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0],
    [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0],
    [-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
]

# The directory containing the data
DIR = "tabledataFinal5Merged/"

tbl = table()

# Utility function to generate a string formatted like the columns
def getColumn(alpha, epsilon, theta, h):
    return f"{alpha}, {epsilon}, {theta}, {h}"

def addTableParameter(table, name, default, delta, values):

    parameter = tableParameter(
        name, 
        0, 
        default, 
        delta, 
        min(values), 
        min(values), 
        max(values), 
        max(values)
    )
    parameter.setTabulatedValues(values)

    table.pushParameter(parameter)

    return table

# ===============================================================
# Primary header unit
# ===============================================================

tbl.setModelName("JhnsnLmPst")
tbl.setModelUnits("counts")
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
# Parameters
# ---------------------------------------------------------------

tbl = addTableParameter(tbl, "SPIN", 0.998, 0.001, spins)
tbl = addTableParameter(tbl, "HEIGHT", 9.0, 0.1, heights)
tbl = addTableParameter(tbl, "INCLINATION", 65.0, 0.1, inclinations)
tbl = addTableParameter(tbl, "ALPHA13", 0.0, 0.1, alphas)
tbl = addTableParameter(tbl, "EPSILON3", 0.0, 0.1, epsilons)

print(tbl)

# ===============================================================
# Reading the data to store in the table model
# ===============================================================

# Reading the CSV
print("Loading CSVs...")

files = {}

for i in os.listdir(DIR):
    if i[-4:] == ".csv":
        print(i)
        files[i[:-4]] = pd.read_csv(DIR + i)

i=0
j=0

# Looping through the parameters in the same order as for generation
print("Parameter loop...")
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:

                    j+=1

                    try:
                        flux = files[str(a)][getColumn(alpha, epsilon, theta, h)].to_numpy()
                    except:
                        flux = np.zeros(1000)
                    
                    if np.count_nonzero(flux) == 0:
                        i+=1

                    # Storing the spectrum and pushing it to the table alongside its 
                    # parameter combination
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

# Saving the file
# tablefile = DIR + "model.FITS"
tablefile = "/home/brad/Documents/SummerInternship/Code/models/TestModel.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

# Printing the number of failed spectra to compare against the logs
print(f"{100*(i/j)}% Failure")
