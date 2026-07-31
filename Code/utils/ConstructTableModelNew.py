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

The format of the csv is as follows:
    - The filename is the spin value of the metric
    - The column headers are in the format "{alpha}, {epsilon}, {theta}, {h}"
    where alpha is the alpha_13 deformation parameter, epsilon is the epsilon_3 
    deformation parameter, theta is the inclination and h is the corona height.
    - The contents of each column are 1000 rows containing the spectrum of the line profile
    generated from the configuration, between 0 and 3 keV for a rest frame energy of 1.0
"""

def addTableParameter(table, name: str, default: float, delta: float, values: list):
    """
    Add a parameter to the table.
    """
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

# ==========================================================================================
# Setup
# ==========================================================================================

# The values each parameter was varied between
# These must be in ascending order for the table model to be assembled correctly
spins, heights, inclinations, alphas, epsilons = [
    [0.0, 0.11088888888888888, 0.22177777777777777, 0.33266666666666667, 0.44355555555555554, 0.5544444444444444, 0.6653333333333333, 0.7762222222222223, 0.8871111111111111, 0.998],
    [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0],
    [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0],
    [-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
]

# The directory containing the data
DIR = "tabledataFinal5Merged/"

# Instantiating the table
tbl = table()

# ==========================================================================================
# Primary header unit
# ==========================================================================================

tbl.setModelName("JhnsnLmPst")
tbl.setModelUnits("counts")
tbl.setisRedshift(False)
tbl.setisAdditive(True)
tbl.setisError(False)
tbl.setisEscale(True)

# ==========================================================================================
# Defining the energies
# ==========================================================================================

tbl.setEnergies(np.linspace(0., 3., 1001))

# ==========================================================================================
# Defining the parameters
# ==========================================================================================

tbl.setNumIntParams(5)
tbl.setNumAddParams(0)

# ------------------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------------------

tbl = addTableParameter(tbl, "SPIN", 0.998, 0.001, spins)
tbl = addTableParameter(tbl, "HEIGHT", 9.0, 0.1, heights)
tbl = addTableParameter(tbl, "INCLINATION", 65.0, 0.1, inclinations)
tbl = addTableParameter(tbl, "ALPHA13", 0.0, 0.1, alphas)
tbl = addTableParameter(tbl, "EPSILON3", 0.0, 0.1, epsilons)

print(tbl)

# ==========================================================================================
# Reading the data to store in the table model
# ==========================================================================================

# Reading the CSV
print("Loading CSVs...")

files = {}

# Getting a dictionary containing all of the data as dataframes
for i in os.listdir(DIR):
    if i[-4:] == ".csv":
        print(i)
        files[i[:-4]] = pd.read_csv(DIR + i)

# Variables for tracking
i = 0
total = len(spins)*len(heights)*len(inclinations)*len(alphas)*len(epsilons)

# ==========================================================================================
# Adding the spectra to the table
# ==========================================================================================

# Looping through the parameters in the same order as for generation
print("Parameter loop...")
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:

                    # Trying to read the data, getting an array of 0 if no data is available
                    try:
                        flux = files[str(a)][f"{alpha}, {epsilon}, {theta}, {h}"].to_numpy()
                    except:
                        flux = np.zeros(1000)
                    
                    # If the spectrum is all 0, incrementing the error counter
                    if np.count_nonzero(flux) == 0:
                        i+=1

                    # Storing the spectrum/combination and adding to the table
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

# ==========================================================================================
# File handling and error reporting
# ==========================================================================================

# Saving the file
tablefile = DIR + "model.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

# Printing the number of failed spectra to compare against the logs
print(f"{100*(i/total)}% Failure")
