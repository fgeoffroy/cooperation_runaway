from __future__ import division
from random import *
import numpy as np
from Utilities import *
from Parameters import *

class DiscreteNorm:

	def __init__(self):
		self.requestVector = np.zeros(NbOfDiscreteRequests)


	def applyNorm(self, offer):
		if offer >= 1:
			return self.requestVector[NbOfDiscreteRequests - 1]	#Cette norme ne repond qu'a des valeurs d'offre entre 0 et 1. Toutes les offres au dela de 1 ont la meme reponse que 1
		else:
			req = self.requestVector[0]
			i = 1
			while offer > (i / NbOfDiscreteRequests):	#Cette norme ne repond qu'a des valeurs d'offre entre 0 et 1
				req = self.requestVector[i]
				i += 1
			return req


	def mutate(self):
		self.requestVector = vMutateDiscreteNorm(self.requestVector)
