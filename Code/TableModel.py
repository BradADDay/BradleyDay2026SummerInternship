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

spins, heights, inclinations, alphas, epsilons = [
    [0.0, 0.2495, 0.499, 0.7485, 0.998],
    [3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0],
    [5.0, 18.333333333333332, 31.666666666666668, 45.0, 58.333333333333336, 71.66666666666667, 85.0],
    [-0.4, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0],
    [-0.4, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
]

tbl = table()

# Utility function to generate a string formatted like the columns
def getColumn(a, h, theta, alpha, epsilon):
    return f"{a}, {h}, {theta}, {alpha}, {epsilon}"

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

spin = tableParameter("SPIN", 0, 0.998, 0.001, -0.998, -0.998, 0.998, 0.998)
spin.setTabulatedValues(spins)

tbl.pushParameter(spin)

# ---------------------------------------------------------------
# Height
# ---------------------------------------------------------------

height = tableParameter("HEIGHT", 0, 10.0, 0.1, 3.0, 3.0, 15.0, 15.0) 
height.setTabulatedValues(heights)

tbl.pushParameter(height)

# ---------------------------------------------------------------
# Inclination
# ---------------------------------------------------------------

inclination = tableParameter("INCLINATION", 0, 60.0, 0.1, 5.0, 5.0, 85.0, 85.0) 
inclination.setTabulatedValues(inclinations)

tbl.pushParameter(inclination)

# ---------------------------------------------------------------
# Alpha13
# ---------------------------------------------------------------

alpha13 = tableParameter("ALPHA13", 0, 0.0, 0.1, -0.4, -0.4, 10.0, 10.0) 
alpha13.setTabulatedValues(alphas)

tbl.pushParameter(alpha13)

# ---------------------------------------------------------------
# Epsilon3
# ---------------------------------------------------------------

epsilon3 = tableParameter("EPSILON3", 0, 0.0, 0.1, -0.4, -0.4, 10.0, 10.0)
epsilon3.setTabulatedValues(epsilons)

tbl.pushParameter(epsilon3)

# ===============================================================
# Reading the data to store in the table model
# ===============================================================

# Reading the CSV
print("Loading CSVs...")
DIR = "tabledata/"
df0 = pd.read_csv(DIR + "0")
df1 = pd.read_csv(DIR + "1")
df2 = pd.read_csv(DIR + "2")
df3 = pd.read_csv(DIR + "3")
df4 = pd.read_csv(DIR + "4")
df5 = pd.read_csv(DIR + "5")

files = [df0, df1, df2, df3, df4, df5]

i=0

# Looping through the parameters in the same order as for generation
print("Parameter loop...")
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:

                    flux=None

                    for file in files:
                        # Pulling the flux from the CSV 
                        try:
                            flux = file[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
                        except:
                            pass
                    
                    if type(flux) == None:
                        flux = np.zeros(1000)
                        i+=1
                    
                    # Storing the spectrum and pushing it to the table alongside its parameter combination
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

# Saving the file
tablefile = "tabledata/model.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

# Printing the number of failed spectra to compare against the logs
print(i)
