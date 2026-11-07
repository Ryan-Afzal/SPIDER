# examplePWRSpentFuel.py
#
# Author:	Ryan Afzal
# Date:		2026-11-06
#
# Example demonstrating the use of SPIDER to plot the activity and radiotoxicity 
# of PWR spent fuel over time. The composition of a standard Westinghouse PWR 
# spent fuel assembly (in 1977) is tabulated, and this data is used to produce 
# a plot similar to the one here, from Wikimedia Commons: 
# https://commons.wikimedia.org/wiki/File:Spent_nuclear_fuel_decay_sievert.jpg
#
# The assembly in question contains 461.4kg of enriched uranium (3.3% U-235).
# This is normalized from per assembly to per metric ton of heavy metal (MTHM).
# 
# Effective radiotoxicity is E(50) for ingestion by an adult ("reference man") member of the public.
# 
# 
# 
# This example uses the following data:
# 
# Nuclide decay data is taken from ICRP Publication 107:
# ICRP, 2008. Nuclear Decay Data for Dosimetric Calculations. ICRP Publication 107. Ann. ICRP 38 (3).
# 	ISBN 978-0-7020-3475-6. https://www.icrp.org/publication.asp?id=ICRP%20Publication%20107.
#
# Spent fuel composition data is from:
# "Chapter 7" in M.-S. Yim, Nuclear Waste Management. Lecture Notes in Energy 83. Springer Nature B.V. 2022.
# 	ISBN:978-94-024-2106-4 doi:10.1007/978-94-024-2106-4_7.
# 
# 
# 
# You will need to download the ICRP Publication 107 data and place the .NDX file in the same directory.
# You will also need to supply your own radiotoxicity conversion factors.
# I recommend those from ICRP Publication 119:
# 
# Effective radiotoxicity conversion [E(50)] factors can be taken from ICRP Publication 119:
# ICRP, 2012. Compendium of Dose Coefficients based on ICRP Publication 60. ICRP Publication 119. Ann. ICRP 41(Suppl.).
# 	ISBN 978-1-4557-5430-4. https://www.icrp.org/publication.asp?id=ICRP%20Publication%20119.

import math
import matplotlib.pyplot as plt
import numpy as np
import sources
import spider


# ===============
# Load decay data
# ===============

print('Loading data...')

# Initialize decay solver, from ICRP-107.
decaySolver = spider.DecaySolver(sources.getNuclidesFromNDX('ICRP-07.NDX'))

# Dictionary of effective radiotoxicity conversion factors (Sv / Bq) per nuclide.

# TODO: Input a dictionary of coefficients here. It should be in the form:
# { 'Nuclide': coefficient }, ex. {'Pu-239': 1E-8 }
coefficientDict = {}



# ============
# Initial data
# ============

print('Initializing...')

# Initial mass of uranium per assembly, in kg.
PWR_M_U_INIT = 461.4

# Final mass of uranium per assembly, in kg.
PWR_M_U_FINAL = 441.2

# Activation products in fresh PWR spent fuel, in MBq per assembly.
PWR_ACTIVATION_PRODUCTS = {
	'H-3': 1.40E7,
	'C-14': 2.54E4,
	'Fe-55': 1.68E6,
	'Ni-59': 7.84E4,
	'Co-58': 1.17E8,
	'Co-60': 1.34E8,
	'Ni-63': 1.13E7,
	'Sn-119m': 1.14E8,
	'Sb-125': 2.91E7
}

# Fission products in fresh PWR spent fuel, in MBq per assembly.
PWR_FISSION_PRODUCTS = {
	'Se-79': 6.99E3,
	'Kr-85': 1.62E8,
	'Sr-90': 1.29E9,
	'Zr-95': 2.63E10,
	'Nb-95': 2.67E10,
	'Tc-99': 2.26E5,
	'Sn-126': 1.32E4,
	'I-129': 5.26E2,
	'I-131': 1.56E10,
	'Xe-133': 3.12E10,
	'Cs-137': 1.77E9,
	'Ce-144': 2.06E10,
	'Eu-154': 1.81E8,
	'Eu-155': 1.19E8
}

# Minor actinides in fresh PWR spent fuel, in MBq per assembly.
PWR_MINOR_ACTINIDES = {
	'Np-237': 5.19E3,
	'Am-241': 1.86E6,
	'Am-243': 2.82E5,
	'Cm-242': 6.45E8,
	'Cm-244': 2.62E7
}

# Plutonium in fresh PWR spent fuel, in MBq per assembly.
PWR_PLUTONIUM = {
	'Pu-238': 3.64E7,
	'Pu-239': 5.18E6,
	'Pu-240': 8.72E6,
	'Pu-241': 2.12E9,
	'Pu-242': 3.07E4
}

# Uranium in fresh PWR spent fuel, in Bq per assembly.
PWR_URANIUM = {
	'U-234': decaySolver.getActivityFromMass( 0.01*1E-2*PWR_M_U_FINAL, 'U-234'),
	'U-235': decaySolver.getActivityFromMass( 0.79*1E-2*PWR_M_U_FINAL, 'U-235'),
	'U-236': decaySolver.getActivityFromMass( 0.40*1E-2*PWR_M_U_FINAL, 'U-236'),
	'U-238': decaySolver.getActivityFromMass(94.29*1E-2*PWR_M_U_FINAL, 'U-238')
}

# Initial uranium in a PWR fuel assembly, in Bq per assembly.
PWR_INIT_URANIUM = {
	'U-235': decaySolver.getActivityFromMass( 3.3*1E-2*PWR_M_U_INIT, 'U-235'),
	'U-238': decaySolver.getActivityFromMass(96.7*1E-2*PWR_M_U_INIT, 'U-238')
}



# =================
# Calculated inputs
# =================

# Mass of natural uranium required per assembly, in kg.
# For this calculation, we assume the following enrichment factors:
# x_nat  = 0.72wt%
# x_prod = 3.3 wt%
# x_tail = 0.30wt%
PWR_M_U_NAT = PWR_M_U_INIT*((0.033 - 0.0030) / (0.0072 - 0.0030))

# Conversion factor between assembly mass and MTHM.
PER_MTHM = 1000 / PWR_M_U_INIT

# Natural uranium required for a PWR fuel assembly, in Bq per assembly.
NAT_URANIUM = {
	'U-235': decaySolver.getActivityFromMass( 0.72*1E-2*PWR_M_U_NAT, 'U-235'),
	'U-238': decaySolver.getActivityFromMass(99.28*1E-2*PWR_M_U_NAT, 'U-238')
}



# ================
# Generate results
# ================

print('Generating data...')

# Array of time points; ranges from 1 to 1E6, with 32 points per order of magnitude.
times = []
for i in range(0, 6):
	for j in range(0, 32):
		times.append(math.pow(10, i) * math.pow(10, j/32))
times.append(1E6)
times = np.array(times)



def getTotalActivityAndRadiotoxicity(nuclidesDict: dict[str, float], t: float) -> tuple[float, float]:
	"""
	Returns a tuple containing the activity and radiotoxicity of a sample of multiple nuclides.

	Args:
		nuclidesDict (float): The sample at time t = 0.
		t (float): The time after t = 0 at which to calculate activity.

	Returns:
		tuple[float, float]: A tuple of the total activity and radiotoxicity in the sample at time t.
	"""
	totalActivity = 0
	totalRadiotoxicity = 0
	for nuclide in nuclidesDict.keys():
		initialAmount = nuclidesDict[nuclide] / decaySolver.getDecayConstant(nuclide)
		a, d = decaySolver.getTotalActivityAtTime(initialAmount, nuclide, t)
		totalActivity += a
		for key in d.keys():
			if (key in coefficientDict):
				totalRadiotoxicity += d[key]*coefficientDict[key]
	return totalActivity, totalRadiotoxicity



# Result arrays; arr[i][0] is activity, arr[i][1] is radiotoxicity.
resultsNU = np.zeros((len(times), 2))# Natural Uranium

resultsUR = np.zeros((len(times), 2))# Uranium
resultsAP = np.zeros((len(times), 2))# Activation Products
resultsFP = np.zeros((len(times), 2))# Fission Products
resultsMA = np.zeros((len(times), 2))# Minor Actinides
resultsPU = np.zeros((len(times), 2))# Plutonium

resultsFAP = np.zeros((len(times), 2))# Fission & Activation Products
resultsMAP = np.zeros((len(times), 2))# Minor Actinides & Plutonium
resultsTOT = np.zeros((len(times), 2))# Total spent fuel

# Partitioning & Transmutation (P&T) of minor actinides and plutonium
results990 = np.zeros((len(times), 2))# P&T of 99.0% of MA, Pu, and Tc-99
results999 = np.zeros((len(times), 2))# P&T of 99.9% of MA, Pu, and Tc-99



# Activities in MBq are multiplied by 1E6; all are multiplied by PER_MTHM.
for i in range(len(times)):
	t = times[i]*spider.YEAR# Physical time

	resultsNU[i] = getTotalActivityAndRadiotoxicity(NAT_URANIUM, t)

	resultsUR[i] = getTotalActivityAndRadiotoxicity(PWR_URANIUM, t)
	resultsAP[i] = getTotalActivityAndRadiotoxicity(PWR_ACTIVATION_PRODUCTS, t)
	resultsFP[i] = getTotalActivityAndRadiotoxicity(PWR_FISSION_PRODUCTS, t)
	resultsMA[i] = getTotalActivityAndRadiotoxicity(PWR_MINOR_ACTINIDES, t)
	resultsPU[i] = getTotalActivityAndRadiotoxicity(PWR_PLUTONIUM, t)

	resultsNU[i] *= PER_MTHM
	resultsUR[i] *= PER_MTHM
	resultsAP[i] *= PER_MTHM * 1E6
	resultsFP[i] *= PER_MTHM * 1E6
	resultsMA[i] *= PER_MTHM * 1E6
	resultsPU[i] *= PER_MTHM * 1E6

	resultsFAP[i] = resultsAP[i] + resultsFP[i]
	resultsMAP[i] = resultsMA[i] + resultsPU[i]
	resultsTOT[i] = resultsUR[i] + resultsAP[i] + resultsFP[i] + resultsMA[i] + resultsPU[i]
	results990[i] = resultsAP[i] + resultsFP[i] + 1E-2*(resultsMA[i] + resultsPU[i])
	results999[i] = resultsAP[i] + resultsFP[i] + 1E-3*(resultsMA[i] + resultsPU[i])

# Transpose the result arrays for easy plotting.
resultsNU = np.transpose(resultsNU)

resultsUR = np.transpose(resultsUR)
resultsAP = np.transpose(resultsAP)
resultsFP = np.transpose(resultsFP)
resultsMA = np.transpose(resultsMA)
resultsPU = np.transpose(resultsPU)

resultsFAP = np.transpose(resultsFAP)
resultsMAP = np.transpose(resultsMAP)
resultsTOT = np.transpose(resultsTOT)
results990 = np.transpose(results990)
results999 = np.transpose(results999)



# ========
# Plotting
# ========

print('Plotting...')

# Colors
YELLOW = '#DDAA33'
RED = '#BB5566'
BLUE = '#004488'
BLACK = '#000000'

def getPlotTemplate():
	fig, ax = plt.subplots(nrows=1, ncols=1)
	ax.grid(visible=True, which='major', axis='both')
	ax.set_xlim(1E2, 1E6)
	ax.set_xlabel('Time after discharge (y)')
	ax.set_xscale('log')
	ax.set_yscale('log')
	return fig, ax

def getActivityPlot():
	fig, ax = getPlotTemplate()
	ax.set_ylabel('Activity (Bq / MTHM)')
	ax.set_ylim(1E11, 1E16)
	return fig, ax

def getRadiotoxicityPlot():
	fig, ax = getPlotTemplate()
	ax.set_ylabel('Radiotoxicity (Sv / MTHM)')
	ax.set_ylim(1E2, 1E8)
	return fig, ax



# Activity vs year (composition)
fig1, ax1 = getActivityPlot()
ax1.plot(times, resultsTOT[0], BLACK, label='Total')
ax1.plot(times, resultsUR[0], BLACK, linestyle=':', label='Uranium and Daughters')
ax1.plot(times, resultsFAP[0], BLUE, label='Fission and Activation Products')
ax1.plot(times, resultsMAP[0], RED, label='Actinides and Daughters')
ax1.plot(times, resultsNU[0], YELLOW, linestyle=':', label='Natural Uranium')
ax1.hlines(resultsNU[0][-1], times[0], times[-1], YELLOW, label='Uranium Ore')
ax1.legend()
fig1.savefig('example_PWRSpentFuel_fig1.svg')

# Activity vs year (reprocessing)
fig2, ax2 = getActivityPlot()
ax2.plot(times, resultsTOT[0], BLACK, label='Total (unprocessed)')
ax2.plot(times, results990[0], RED, label='P&T of 99.0% of MA + Pu')
ax2.plot(times, results999[0], BLUE, label='P&T of 99.9% of MA + Pu')
ax2.vlines([633.106, 2157.4, 194704], 1E2, resultsNU[0][-1], colors=[BLUE, RED, BLACK], linestyles='--')
ax2.hlines(resultsNU[0][-1], times[0], times[-1], YELLOW, label='Uranium Ore')
ax2.legend()
fig2.savefig('example_PWRSpentFuel_fig2.svg')

fig3, ax3 = getRadiotoxicityPlot()
ax3.plot(times, resultsTOT[1], BLACK, label='Total')
ax3.plot(times, resultsUR[1], BLACK, linestyle=':', label='Uranium and Daughters')
ax3.plot(times, resultsFAP[1], BLUE, label='Fission and Activation Products')
ax3.plot(times, resultsMAP[1], RED, label='Actinides and Daughters')
ax3.plot(times, resultsNU[1], YELLOW, linestyle=':', label='Natural Uranium')
ax3.hlines(resultsNU[1][-1], times[0], times[-1], YELLOW, label='Uranium Ore')
ax3.legend()
fig3.savefig('example_PWRSpentFuel_fig3.svg')

fig4, ax4 = getRadiotoxicityPlot()
ax4.plot(times, resultsTOT[1], BLACK, label='Total (unprocessed)')
ax4.plot(times, results990[1], RED, label='P&T of 99.0% of MA + Pu')
ax4.plot(times, results999[1], BLUE, label='P&T of 99.9% of MA + Pu')
ax4.vlines([281.062, 562.720, 115182], 1E2, resultsNU[1][-1], colors=[BLUE, RED, BLACK], linestyles='--')
ax4.hlines(resultsNU[1][-1], times[0], times[-1], YELLOW, label='Uranium Ore')
ax4.legend()
fig4.savefig('example_PWRSpentFuel_fig4.svg')



# Replication of the plot from Wikimedia Commons
fig5, ax5 = getRadiotoxicityPlot()
ax5.plot(times, resultsTOT[1], 'black', label='Total')
ax5.plot(times, resultsUR[1], 'black', label='Uranium')
ax5.plot(times, resultsFAP[1], 'purple', label='Fission & Activation Products')
ax5.plot(times, resultsPU[1], 'blue', label='Plutonium')
ax5.plot(times, resultsMA[1], 'red', label='Minor Actinides')
ax5.plot(times, resultsNU[1], 'k:', label='Natural U')
ax5.hlines(resultsNU[1][-1], 1, 1E6, 'green', label='Uranium Ore')
ax5.legend()
fig5.savefig('example_PWRSpentFuel_fig5.svg')



# Uncomment to show the plots at the end.
#plt.show()

print('Finished.')
