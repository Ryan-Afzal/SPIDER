# spider.py
#
# Author:	Ryan Afzal
# Date:		2026-10-17
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
		"""Returns the half-life of nuclide n, in s."""
		return self._nuclidesDict[n][2]

	def getDecayConstant(self, n: str) -> float:
		"""Returns the decay constant of nuclide n, in decay/s."""
		return decayConst(self.getHalfLife(n))

	def getSpecificActivity(self, n: str) -> float:
		"""Returns the specific activity of nuclide n, in Bq/kg."""
		return self.getDecayConstant(n) * NA / self._nuclidesDict[n][1] * 1000

	def getAtomsFromMass(self, mass: float, n: str) -> float:
		"""Returns the number of atoms of a mass, in kg, of nuclide n."""
		return NA / self._nuclidesDict[n][1] * mass * 1000

	def getActivityFromMass(self, mass: float, n: str) -> float:
		"""Returns the activity of a mass of nuclide n, in Bq."""
		return self.getDecayConstant(n) * self.getAtomsFromMass(mass, n)

	def _getMatrix(self, constants: list[float], fractions: list[float]) -> list[list[float]]:
		"""Get the coefficient matrix for the provided nuclei decay constant chain and branching fractions."""

		if (len(constants) != (len(fractions) + 1)):
			raise ValueError(f'Constants list must contain 1 more element than the fractions list ({len(constants)}, {len(fractions)})')
		
		matrix = []

		num = 1
		for i in range(len(constants) - 1):
			num *= (constants[i] * fractions[i])

		for i in range(len(constants)):
			den = 1
			for j in range(len(constants)):
				if (i != j):
					den *= constants[j] - constants[i]
			c = num / den
			matrix.append([constants[i], c])

		return matrix

	def getTotalActivityAtTime(self, N0: float, nuclide: str, t: float) -> tuple[float, dict[str, float]]:
		"""
		Returns the total activity, at time t, of the specified starting nuclide activity.

		Args:
			N0 (float): The initial quantity of material in atoms.
			nuclide (str): The starting nuclide of activity A0. E.g. 'Pu-239'.
			t (float): The desired time at which to calculate the activity.
		Returns:
			tuple[float, dict[str, float]]: A tuple containing the total activity, and a dict containing the activities of each daughter nuclide.
		"""
		totalActivities = {}
		totalActivity = self._getTotalActivityAtTime(N0, [], [], nuclide, totalActivities, t)

		return totalActivity, totalActivities

	def _getTotalActivityAtTime(self, P0: float, decayChain: list[float], branchFractions: list[float], nuclide: str, activities: dict[str,float], t: float) -> float:
		if (nuclide == 'SF'):
			# TODO
			return 0
		elif (nuclide not in self._nuclidesDict):
			# Stable nuclide, so no activity
			return 0

		nuclideDecayConst = self.getDecayConstant(nuclide)
		decayChain.append(nuclideDecayConst)

		matrix = self._getMatrix(decayChain, branchFractions)
		activity = 0
		for k in range(len(matrix)):
			activity += matrix[k][1]*math.exp(-matrix[k][0]*t)
		activity *= (nuclideDecayConst*P0)

		if (nuclide in activities):
			activities[nuclide] += activity
		else:
			activities[nuclide] = activity

		for daughter in self._nuclidesDict[nuclide][3]:
			newDecayChain = decayChain.copy()
			newBranchFractions = branchFractions.copy()
			newBranchFractions.append(daughter[1])
			activity += self._getTotalActivityAtTime(P0, newDecayChain, newBranchFractions, daughter[0], activities, t)

		return activity
