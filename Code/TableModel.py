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

spins = [0.0, 0.16633333333333333, 0.33266666666666667, 0.499, 0.6653333333333333, 0.8316666666666667, 0.998]
spin = tableParameter("SPIN", 0, 0.998, 0.001, -0.998, -0.998, 0.998, 0.998)
spin.setTabulatedValues(spins)

tbl.pushParameter(spin)

# ---------------------------------------------------------------
# Height
# ---------------------------------------------------------------

heights = [3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0]
height = tableParameter("HEIGHT", 0, 10.0, 0.1, 3.0, 3.0, 15.0, 15.0) 
height.setTabulatedValues(heights)

tbl.pushParameter(height)

# ---------------------------------------------------------------
# Inclination
# ---------------------------------------------------------------

inclinations = [5.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0]
inclination = tableParameter("INCLINATION", 0, 60.0, 0.1, 5.0, 5.0, 85.0, 85.0) 
inclination.setTabulatedValues(inclinations)

tbl.pushParameter(inclination)

# ---------------------------------------------------------------
# Alpha13
# ---------------------------------------------------------------

alphas = [-0.4, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
alpha13 = tableParameter("ALPHA13", 0, 0.0, 0.1, -0.4, -0.4, 10.0, 10.0) 
alpha13.setTabulatedValues(alphas)

tbl.pushParameter(alpha13)

# ---------------------------------------------------------------
# Epsilon3
# ---------------------------------------------------------------

epsilons = [-0.4, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
epsilon3 = tableParameter("EPSILON3", 0, 0.0, 0.1, -0.4, -0.4, 10.0, 10.0)
epsilon3.setTabulatedValues(epsilons)

tbl.pushParameter(epsilon3)

# ===============================================================
# Reading the data to store in the table model
# ===============================================================

# Reading the CSV
print("Loading CSVs...")
df2 = pd.read_csv("data/spectra2.csv")
df3 = pd.read_csv("data/spectra3.csv")
df4 = pd.read_csv("data/NegativeDeformation.csv")
i=0

# Looping through the parameters in the same order as for generation
print("Parameter loop...")
for a in spins:
    for h in heights:
        for theta in inclinations:
            for alpha in alphas:
                for epsilon in epsilons:

                    # Pulling the flux from the CSV 
                    try:
                        flux = df2[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
                    except:
                        try:
                            flux = df3[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
                        except:
                            try:
                                flux = df4[getColumn(a, h, theta, alpha, epsilon)].to_numpy()
                            # Setting the flux to zero if the spectrum failed
                            # This may be better done through interpolation in the future, 
                            # however it was often consecutive parameter combinations that failed
                            except:
                                flux = np.zeros(1000)
                    
                    # Storing the spectrum and pushing it to the table alongside its parameter combination
                    spec = tableSpectrum()
                    spec.setParameterValues(np.array([a, h, theta, alpha, epsilon]))
                    spec.setFlux(flux)
                    tbl.pushSpectrum(spec)

# Saving the file
tablefile = "model2.FITS"
if (os.path.exists(tablefile)): 
    os.remove(tablefile)
status = tbl.write(tablefile)
if status != 0: 
    print("Failed to write test.mod: status = ", status)

# Printing the number of failed spectra to compare against the logs
print(i)
