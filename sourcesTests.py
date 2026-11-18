# sourcesTests.py
#
# Author:	Ryan Afzal
# Date:		2026-10-15
#
# Unit test script to verify sources.py using the ICRP Publication 107 dataset.
#
# To test, place the ICRP-07.NDX file in the same directory and run.

import src.sources as sources
import unittest

class TestNDX(unittest.TestCase):
	_nuclidesDict = {}
	def setUp(self):
		if (len(self._nuclidesDict) == 0):
			self._nuclidesDict = sources.getNuclidesFromNDX('ICRP-07.NDX')
	def test_loadndx(self):
		"""Tests if anything has been loaded."""
		self.assertGreater(len(self._nuclidesDict), 0)
	def test_loadICRP07_numNuclides(self):
		"""Tests if all nuclides (no more, no fewer) have been loaded."""
		self.assertEqual(len(self._nuclidesDict), 1252)
	def test_loadICRP07_halfLifeSingle(self):
		"""Tests if U-233 has the correct half-life loaded, with its single-letter unit suffix handled properly."""
		self.assertEqual(self._nuclidesDict['U-233'][2], 5020531200000)
	def test_loadICRP07_halfLifeDouble(self):
		"""Tests if Ra-219 has the correct half-life loaded, with its double-letter unit suffix handled properly."""
		self.assertEqual(self._nuclidesDict['Ra-219'][2], 10E-3)
	def test_loadICRP07_decayModes(self):
		"""Tests if Cf-246 has been loaded correctly with the correct atomic mass, which implies decay modes are handled correctly."""
		self.assertEqual(self._nuclidesDict['Cf-246'][1], 246.068805)
	def test_loadICRP07_threeDaughters(self):
		"""Tests if Ac-226 is loaded properly, with its three daughter nuclides loaded with correct branching fractions."""
		daughters = self._nuclidesDict['Ac-226'][3]
		self.assertEqual(len(daughters), 3)
		self.assertEqual(daughters[0][0], 'Th-226')
		self.assertEqual(daughters[0][1], 8.3000E-01)
		self.assertEqual(daughters[1][0], 'Ra-226')
		self.assertEqual(daughters[1][1], 1.7000E-01)
		self.assertEqual(daughters[2][0], 'Fr-222')
		self.assertEqual(daughters[2][1], 6.0000E-05)
	def test_loadICRP07_fourDaughters(self):
		"""Tests if Es-254m is loaded properly, with its four* daughter nuclides loaded with correct branching fractions (spontaneous fission is counted as a daughter)."""
		daughters = self._nuclidesDict['Es-254m'][3]
		self.assertEqual(len(daughters), 4)
		self.assertEqual(daughters[0][0], 'Fm-254')
		self.assertEqual(daughters[0][1], 9.8000E-01)
		self.assertEqual(daughters[1][0], 'Bk-250')
		self.assertEqual(daughters[1][1], 3.2000E-03)
		self.assertEqual(daughters[2][0], 'Cf-254')
		self.assertEqual(daughters[2][1], 7.6000E-04)
		self.assertEqual(daughters[3][0], 'SF')
		self.assertEqual(daughters[3][1], 4.5000E-04)

if __name__ == "__main__":
	unittest.main()