# spider.py
#
# Author:	Ryan Afzal
# Date:		2026-10-15
#
# Contains various methods to generate and analyze nuclide decay chains.

import math

NA = 6.022E23
CURIE = 3.7E10
HOUR = 3600
DAY = 24*HOUR
YEAR = 365*DAY

def decayConst(t12):
	return math.log(2)/t12

class DecaySolver:

	def __init__(self, nuclidesDict = {}):
		self._nuclidesDict = nuclidesDict

	def getHalfLife(self, n: str) -> float:
		return self._nuclidesDict[n][2]

	def getDecayConstant(self, n: str) -> float:
		return decayConst(self.getHalfLife(n))

	def getSpecificActivity(self, n: str) -> float:
		"""Returns the specific activity of nuclide n in Bq/kg."""
		return self.getDecayConstant(n) * NA / self._nuclidesDict[n][1] * 1000

	def getAtomsFromMass(self, mass: float, n: str) -> float:
		"""Returns the number of atoms of a mass in kg of nuclide."""
		return NA / self._nuclidesDict[n][1] * mass * 1000

	def getActivityFromMass(self, mass: float, n: str) -> float:
		"""Returns the activity of a mass of nuclide n in Bq."""
		return self.getDecayConstant(n) * self.getAtomsFromMass(mass, n)

	