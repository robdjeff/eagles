# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 12:41:08 2022

@author: Rob Jeffries

Produces and plots a set of isochrones of EWLi vs Teff (and the modelled dispersion)
for a defined list of ages; and saves them as individual text files.
Optionally it can plot the model dispersion around each isochrone.

"""

import numpy as np
import matplotlib.pyplot as plt
# import the EWLi prediction model from the main EAGLES code
from eagles import AT2EWm
from eagles import eAT2EWm


# Modify the list below for the isochrones you want to produce (units are Myr)
ages = [5, 10, 15, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000]

# e.g. for fewer isochrones
#ages = [10, 100, 1000]

# Modify this flag to true if you want the plot to include the dispersion
# look smessy if there are many closely spaced isochrones
plot_dispersion = False

# The step in logarithmic temperature
tstep = 0.002

# parameters for the plot
plt.xlabel('Teff (K)')
plt.ylabel('LiEW (mA)')
plt.xlim(6500, 3000)

# set up a an equally spaced set of log temperatures between 3000 and 6500 K
lteff = np.arange(3.4772, 3.8130, tstep)

# loop over the ages
for t in ages :

    lAge = np.log10(t)+6  # log age in years
    ewm = AT2EWm(lteff, lAge)
    eewm = eAT2EWm(lteff, lAge)

    # save the results as a simple .txt file    
    name = 'iso_'+str(t)+'.txt'
    np.savetxt(name, np.column_stack((10**lteff, ewm, eewm)), fmt='%.1f %.1f %.1f', delimiter=' ', header = "Teff(K) EWLim(mA) eEWLi(mA)")
 
    plt.plot(10**lteff, ewm)

    # if the plot_dispersion flag then shade the dispersion region
    # looks quite messy if there are lots of isochrones because of the overlap
    if plot_dispersion :
        plt.fill_between(10**lteff, ewm-eewm, ewm+eewm, alpha=0.3)

# For the default list of ages and plot_dispersion = False, this is Fig.2 from the paper    
plt.show()

