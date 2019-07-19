from __future__ import division
from random import *
import numpy as np
from math import *
from Utilities import *
from Parameters import *

class PolynomialNorm:

	def __init__(self):
		self.polynomeCoefficients = np.zeros(NbOfPolynomialCoefficients)


	def applyNorm(self, offer):
		req = 0
		i = 0
		while i < NbOfPolynomialCoefficients:
			req += self.polynomeCoefficients[i] * pow(offer, i)
			i += 1
		return req


	def mutate(self):
		self.polynomeCoefficients = vMutatePolynomialNorm(self.polynomeCoefficients)
