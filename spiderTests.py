# spiderTests.py
#
# Author:	Ryan Afzal
# Date:		2026-10-17
#
# Unit test script to verify spider.py

import spider

import math
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
		'Pa-234m': ['Pa-234m', 234.043306, 1.159*60, [['U-234', 0.9984], ['Pa-234', 0.0016]]],
		'U-234': ['U-234', 234.0409503, 245500*365*24*3600, [['Th-230', 1]]],
		'Pa-234': ['Pa-234', 234.0433056, 6.70*3600, [['U-234', 1]]]
	}
	
	def test_halfLife(self):
		"""Tests getHalfLife()."""
		decaySolver = spider.DecaySolver(self._nuclides)
		self.assertEqual(decaySolver.getHalfLife('Pa-234m'), 1.159*60)
	def test_specificActivity(self):
		"""Tests getSpecificActivity() to within 5 decimal places."""
		decaySolver = spider.DecaySolver(self._nuclides)
		sa = math.log(2)/(4.463E9*365*24*3600) * 1/238.050787 * 1000 * spider.NA
		self.assertEqual(decaySolver.getSpecificActivity('U-238'), sa)
	def test_getAtomsFromMass(self):
		"""Tests getAtomsFromMass()."""
		decaySolver = spider.DecaySolver(self._nuclides)
		m = 6.7
		a = m * 1000 / 238.050787 * 6.022E23
		self.assertEqual(decaySolver.getAtomsFromMass(m, 'U-238'), a)
	def test__getMatrix(self):
		"""
		Tests _getMatrix() using the U-238 branching through Pa-234.
		It calculates the matrix for U-238->Th-234->Pa-234m->Pa-234->U-234,
		which would be used to determine the activity of U-234.
		"""
		decaySolver = spider.DecaySolver(self._nuclides)

		# Decay chain
		c = [
			decaySolver.getDecayConstant('U-238'),
			decaySolver.getDecayConstant('Th-234'),
			decaySolver.getDecayConstant('Pa-234m'),
			decaySolver.getDecayConstant('Pa-234'),
			decaySolver.getDecayConstant('U-234')
		]

		# Branching fractions (fraction going to the next element)
		b = [
			1,
			1,
			0.0016,
			1
		]

		num = c[0]*b[0]*c[1]*b[1]*c[2]*b[2]*c[3]*b[3]
		expectedMatrix = [
			[c[0], 
			 num / (c[1] - c[0]) / (c[2] - c[0]) / (c[3] - c[0]) / (c[4] - c[0])
			 ],
			[c[1],
			 num / (c[0] - c[1]) / (c[2] - c[1]) / (c[3] - c[1]) / (c[4] - c[1])
			 ],
			[c[2],
			 num / (c[0] - c[2]) / (c[1] - c[2]) / (c[3] - c[2]) / (c[4] - c[2])
			 ],
			[c[3],
			 num / (c[0] - c[3]) / (c[1] - c[3]) / (c[2] - c[3]) / (c[4] - c[3])
			 ],
			[c[4],
			 num / (c[0] - c[4]) / (c[1] - c[4]) / (c[2] - c[4]) / (c[3] - c[4])
			 ]
		]

		matrix = decaySolver._getMatrix(c, b)

		for i in range(len(matrix)):
			for j in range(2):
				self.assertAlmostEqual(matrix[i][j], expectedMatrix[i][j], 7)
	def test_getTotalActivityAtTime1(self):
		"""Tests getTotalActivityAtTime() using Pa-234 at various times."""
		decaySolver = spider.DecaySolver(self._nuclides)
		a0 = spider.CURIE

		l0 = decaySolver.getDecayConstant('Pa-234')
		l1 = decaySolver.getDecayConstant('U-234')
		n0 = a0/l0

		def predActivityPa234(t):
			return a0*math.exp(-l0*t)

		def predActivityU234(t):
			c0 = l0 / (l1 - l0)
			c1 = l0 / (l0 - l1)
			return l1*n0*(c0*math.exp(-l0*t) + c1*math.exp(-l1*t))

		times = [1E-3, 1E-2, 1E-1, 1, 1E1, 1E2, 1E3, 1E4, 1E5]
		for time in times:
			t = time*spider.YEAR
			total, d = decaySolver.getTotalActivityAtTime(n0, 'Pa-234', t)
			self.assertEqual(d['Pa-234'], predActivityPa234(t))
			self.assertEqual(d['U-234'], predActivityU234(t))
	def test_getTotalActivityAtTime2(self):
		"""Tests getTotalActivityAtTime() using Pa-234 from U-238 at various points in time."""
		decaySolver = spider.DecaySolver(self._nuclides)
		a0 = spider.CURIE

		# Decay chain
		c = [
			decaySolver.getDecayConstant('U-238'),
			decaySolver.getDecayConstant('Th-234'),
			decaySolver.getDecayConstant('Pa-234m'),
			decaySolver.getDecayConstant('Pa-234')
		]
		
		# Branching fractions (fraction going to the next element)
		b = [
			1,
			1,
			0.0016
		]

		c1 = [
			c[0],
			c[1],
			c[2],
			c[3],
			decaySolver.getDecayConstant('U-234')
		]
		b1 = [
			1,
			1,
			0.0016,
			1
		]
		c2 = [
			c[0],
			c[1],
			c[2],
			c1[4]
		]
		b2 = [
			1,
			1,
			0.9984
		]

		n0 = a0 / c[0]

		matrix = decaySolver._getMatrix(c, b)
		matrix1 = decaySolver._getMatrix(c1, b1)
		matrix2 = decaySolver._getMatrix(c2, b2)

		def predActivityPa234(t):
			activity = 0
			for item in matrix:
				activity += item[1]*math.exp(-item[0]*t)
			activity *= (c[-1]*n0)
			return activity

		def predActivityU234(t):
			return predActivityU2341(t) + predActivityU2342(t)
		def predActivityU2341(t):
					activity = 0.0
					for item in matrix1:
						activity += item[1]*math.exp(-item[0]*t)
					activity *= (c1[-1]*n0)
					return activity
		def predActivityU2342(t):
					activity = 0.0
					for item in matrix2:
						activity += item[1]*math.exp(-item[0]*t)
					activity *= (c2[-1]*n0)
					return activity

		times = [1E-3, 1E-2, 1E-1, 1, 1E1, 1E2, 1E3, 1E4, 1E5]
		for time in times:
			t = time*spider.YEAR
			total, d = decaySolver.getTotalActivityAtTime(n0, 'U-238', t)
			self.assertEqual(d['Pa-234'], predActivityPa234(t))
			self.assertEqual(d['U-234'], predActivityU234(t))

if __name__ == "__main__":
	unittest.main()