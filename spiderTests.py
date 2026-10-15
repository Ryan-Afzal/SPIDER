# spiderTests.py
#
# Author:	Ryan Afzal
# Date:		2026-10-15
#
# Unit test script to verify spider.py

import spider

import unittest

class TestSPIDER(unittest.TestCase):
	def test_decayConst(self):
		"""Tests that the decayConst() method works properly on Pu-239. Tests output to within 5 decimal places."""
		t12 = 24110*365*24*3600
		self.assertAlmostEqual(spider.decayConst(t12), 9.11636371e-13, 5)

class TestDecaySolver(unittest.TestCase):
	_nuclides = {
		'H-3': ['H-3', 3.01604928, 12.32*365*24*3600, [['He-3', 1]]],
		'U-238': ['U-238', 238.050787, 4.463E9*365*24*3600, [['Th-234', 1]]],
		'Th-234': ['Th-234', 234.043600, 24.107*3600*24, [['Pa-234m', 1]]],
		'Pa-234m': ['Pa-234m', 234.043306, 1.159*60, [['U-234', 0.9984], ['Pa-234', 0.0016]]]
	}
	
	def test_halfLife(self):
		"""Tests getHalfLife()."""
		decaySolver = spider.DecaySolver(self._nuclides)
		self.assertEqual(decaySolver.getHalfLife('Pa-234m'), 1.159*60)
	def test_specificActivity(self):
		"""Tests getSpecificActivity() to within 5 decimal places."""
		decaySolver = spider.DecaySolver(self._nuclides)
		sa = 1/238.050787 * 1000 * spider.NA * 0.69314718056/(4.463E9*3600*24*365)
		self.assertAlmostEqual(decaySolver.getSpecificActivity('U-238'), sa, 5)
	def test_getAtomsFromMass(self):
		"""Tests getAtomsFromMass()."""
		decaySolver = spider.DecaySolver(self._nuclides)
		m = 6.7
		a = m * 1000 / 238.050787 * 6.022E23
		self.assertEqual(decaySolver.getAtomsFromMass(m, 'U-238'), a)
	

if __name__ == "__main__":
	unittest.main()