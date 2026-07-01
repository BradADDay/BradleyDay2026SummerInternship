from heasp import *
import numpy as np
import os
import pandas as pd

<<<<<<< HEAD
"""
Convert the csv generated from DataTable.jl into a FITS file following the XSPEC
table model specifications. 

See https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/ogip_92_009/ogip_92_009.pdf
for more information.

This makes use of the python wrapper for the heasp C module.

See https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/headas/heasp/node1.html
for more information.
"""

tbl = table()

# Utility function to generate a string formatted like the columns
def getColumn(a, h, theta, alpha, epsilon):
    return f"{a}, {h}, {theta}, {alpha}, {epsilon}"

# ===============================================================
# Primary header unit
# ===============================================================
=======
tbl = table()

def getColumn(a, h, theta, alpha, epsilon):
    return f"{a}, {h}, {theta}, {alpha}, {epsilon}"

# Primary header
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125

tbl.setModelName("JhnsnLmPst")
tbl.setModelUnits(" ")
tbl.setisRedshift(False)
tbl.setisAdditive(True)
tbl.setisError(False)
tbl.setisEscale(True)

<<<<<<< HEAD
# ===============================================================
# Defining the energies
# ===============================================================

tbl.setEnergies(np.linspace(0.2, 1.5, 1001))

# ===============================================================
# Defining the parameters
# ===============================================================
=======
# Energies

tbl.setEnergies(np.linspace(0.2, 1.5, 1001))

# Parameters
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125

tbl.setNumIntParams(5)
tbl.setNumAddParams(0)

<<<<<<< HEAD
# ---------------------------------------------------------------
# Spin
# ---------------------------------------------------------------

=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
spins = [-0.998, -0.6653333333333333, -0.33266666666666667, 0.0, 0.33266666666666667, 0.6653333333333333, 0.998]
spin = tableParameter("SPIN", 0, 0.998, 0.001, -0.998, -0.998, 0.998, 0.998)
spin.setTabulatedValues(spins)

tbl.pushParameter(spin)

<<<<<<< HEAD
# ---------------------------------------------------------------
# Height
# ---------------------------------------------------------------

=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
heights = [1.5, 4.666666666666667, 7.833333333333333, 11.0, 14.166666666666666, 17.333333333333332, 20.5, 23.666666666666668, 26.833333333333332, 30.0]
height = tableParameter("HEIGHT", 0, 10.0, 0.1, 1.5, 1.5, 30.0, 30.0) 
height.setTabulatedValues(heights)

tbl.pushParameter(height)

<<<<<<< HEAD
# ---------------------------------------------------------------
# Inclination
# ---------------------------------------------------------------

=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
inclinations = [5.0, 13.88888888888889, 22.77777777777778, 31.666666666666668, 40.55555555555556, 49.44444444444444, 58.333333333333336, 67.22222222222223, 76.11111111111111, 85.0]
inclination = tableParameter("INCLINATION", 0, 60.0, 0.1, 5.0, 5.0, 85.0, 85.0) 
inclination.setTabulatedValues(inclinations)

tbl.pushParameter(inclination)

<<<<<<< HEAD
# ---------------------------------------------------------------
# Alpha13
# ---------------------------------------------------------------

=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
alphas = [0.0, 8.333333333333334, 16.666666666666668, 25.0, 33.333333333333336, 41.666666666666664, 50.0]
alpha13 = tableParameter("ALPHA13", 0, 0.0, 0.1, 0.0, 0.0, 50.0, 50.0) 
alpha13.setTabulatedValues(alphas)

tbl.pushParameter(alpha13)

<<<<<<< HEAD
# ---------------------------------------------------------------
# Epsilon3
# ---------------------------------------------------------------

=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
epsilons = [0.0, 3.3333333333333335, 6.666666666666667, 10.0, 13.333333333333334, 16.666666666666668, 20.0, 23.333333333333332, 26.666666666666668, 30.0]
epsilon3 = tableParameter("EPSILON3", 0, 0.0, 0.1, 0.0, 0.0, 30.0, 30.0)
epsilon3.setTabulatedValues(epsilons)

tbl.pushParameter(epsilon3)

<<<<<<< HEAD
# ===============================================================
# Reading the data to store in the table model
# ===============================================================

# Reading the CSV
df = pd.read_csv("data/spectra.csv")
i=0

# Looping through the parameters in the same order as for generation
=======
df = pd.read_csv("data/spectra.csv")

i=0

>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:
<<<<<<< HEAD

                    # Pulling the flux from the CSV 
                    try:
                        flux = df[getColumn(a, h, theta, alpha, epsilon)].to_numpy()

                    # Setting the flux to zero if the spectrum failed
                    # This may be better done through interpolation in the future, 
                    # however it was often consecutive parameter combinations that failed
=======
                    try:
                        flux = df[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
                    except:
                        flux = np.zeros(1000)
                        i+=1
                    
<<<<<<< HEAD
                    # Storing the spectrum and pushing it to the table alongside its parameter combination
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

<<<<<<< HEAD
# Saving the file
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
tablefile = "test.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

<<<<<<< HEAD
# Printing the number of failed spectra to compare against the logs
=======
>>>>>>> c482177c9ab188292c53cf672ee3e0f595a34125
print(i)