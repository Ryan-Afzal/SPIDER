# sources.py
#
# Author:	Ryan Afzal
# Date:		2026-10-15
#
# Various methods for retreiving a nuclide dictionary.

from typing import TypedDict

def getNuclidesFromNDX(path: str) -> dict:
	"""
	Pulls nuclide data from the International Commission on Radiological Protection (ICRP)'s Publication 107 dataset.
	The data is in the form of a .NDX file, "ICRP-07.NDX".

	Data is sourced from:
	ICRP, 2008. Nuclear Decay Data for Dosimetric Calculations. ICRP Publication 107. Ann. ICRP 38 (3). ISBN 978-0-7020-3475-6. https://www.icrp.org/publication.asp?id=ICRP%20Publication%20107.

	Args:
		path (str): Path to the provided .NDX file.
	
	Returns:
		dict: Dictionary of nuclides, indexed by name, containing their name, atomic mass, half-life, and list of daughter isotopes.
	"""
	nuclidesDict = {}

	# NDX files use the following format:
	# 	Nuclide (parent) name; ex. Pu-239
	#	Half-life with unit; ex. 2.39E+4y
	#	(unused) Decay modes; A (alpha), B- (beta-), B+ (beta+), EC (electron capture), IT (internal transition), SF (spontaneous fission)
	#	(unused) Pointer-1 I7 Location of nuclide in ICRP-07.RAD file
	#	(unused) Pointer-2 I7 Location of nuclide in ICRP-07.BET file
	#	(unused) Pointer-3 I7 Location of nuclide in ICRP-07.ACK file
	#	(unused) Pointer-4 I6 Location of nuclide in ICRP-07.NSF file
	#	A block of fields, repeated for each daughter nuclide:
	#	(unused) Energy of emitted alphas in MeV/transition
	#	(unused) Energy of emitted electrons in MeV/nuclear transition
	#	(unused) Energy of emitted photons in MeV/nuclear transition
	#	(unused) Number of photons of energy less than 10 keV per nuclear transition
	#	(unused) Number of photons of energy greater than 10 keV per nuclear transition
	#	(unused) Number of beta transitions per nuclear transition
	#	(unused) Number of monoenergetic electrons per nuclear transition
	#	(unused) Number of alpha transitions per nuclear transition
	#	Atomc mass in amu
	#	(unused) Air kerma-rate constant (Gy-m2 / Bq-s)
	#	(unused) Air kerma coefficient (Gy-m2 / Bq-s)

	with open(path, 'r') as inputFile:
		# NDX files have a header line, so we skip this.
		inputFile.readline()

		# Loop through each line in the index and pull out isotope decay data.
		for line in inputFile.readlines():
			# Blot out the "decay modes" section, because it sometimes meets its maximum width and this screws up split().
			line = line[0:17] + '########' + line[25:]

			# Separate the half-life from the decay modes.
			hlUnitIndex = line.find('ms')
			if (hlUnitIndex != -1):
				line = line[0:hlUnitIndex+2] + ' ' + line[hlUnitIndex+2:]
			else:
				hlUnitIndex = line.find('us')
				if (hlUnitIndex != -1):
					line = line[0:hlUnitIndex+2] + ' ' + line[hlUnitIndex+2:]
				else:
					hlUnitIndex = line.find('µs')
					if (hlUnitIndex != -1):
						line = line[0:hlUnitIndex+2] + ' ' + line[hlUnitIndex+2:]

			# Split values by whitespace.
			values = line.split()

			# Nuclide name (ex. Pu-239).
			nuclideName = values[0]

			# Half-life has an extra character or two for the unit.
			if (values[1][-2].isdigit()):
				nuclideHalfLife = float(values[1][0:-1])
				halfLifeUnit = values[1][-1]
			else:
				nuclideHalfLife = float(values[1][0:-2])
				halfLifeUnit = values[1][-2:]
			halfLifeMult = -1
			match halfLifeUnit:
				case 'µs':
					halfLifeMult = 1E-6
				case 'us':
					halfLifeMult = 1E-6
				case 'ms':
					halfLifeMult = 1E-3
				case 's':
					halfLifeMult = 1
				case 'm':
					halfLifeMult = 60
				case 'h':
					halfLifeMult = 3600
				case 'd':
					halfLifeMult = 3600*24
				case 'y':
					halfLifeMult = 3600*24*365
			if (halfLifeMult == -1):
				raise ValueError(f'Half-life unit \"{halfLifeUnit}\" for nuclide {nuclideName} is invalid.')
			# Half-life in seconds
			nuclideHalfLife *= halfLifeMult

			# The last two values are sometimes formatted without a space between them.
			if (len(values[-1]) == 18):
				swapVal = values[-1][-9:]
				values[-1] = values[-1][0:-9]
				values.append(swapVal)
			elif (len(values[-1]) > 18):
				raise ValueError(f'Malformed final value \"{values[-1]}\" for nuclide {nuclideName}.')

			# Atomic mass in amu (Da).
			nuclideAtomicMass = float(values[-3])

			nuclideDaughters = []
			index = 7
			while (values[index] != '0' and index < 17):
				branchFraction = float(values[index + 2])
				nuclideDaughters.append([values[index], branchFraction])
				index += 3

			# Add this line's entry to the dictionary.
			nuclidesDict[nuclideName] = [nuclideName, nuclideAtomicMass, nuclideHalfLife, nuclideDaughters]

	# Return the final dictionary.
	return nuclidesDict