from __future__ import division
from math import *
from random import *
from Utilities import *
from Parameters import *

if IsReactionNorm:
	if IsDiscreteNorm:
		from DiscreteNorm import *
	elif IsPolynomialNorm:
		from PolynomialNorm import *

class Individual:

	def __init__(self):
		self._offerGene = InitialOffer
		self._offer = InitialOffer
		if IsReactionNorm:
			if IsDiscreteNorm:
				self.reactionNorm = DiscreteNorm()
			elif IsPolynomialNorm:
				self.reactionNorm = PolynomialNorm()
			self._request = self.reactionNorm.applyNorm(self.offer)
		else:
			self._request = InitialRequest
		self._payoff = 0
		self._cumulatedPayoff = 0
		self._fecundity = 0
		self._partnerIndex = -1


	def reset(self):
		self.partnerIndex = -1
		self.payoff = 0
		self.cumulatedPayoff = 0
		#Apply phenotypic noise
		if NoiseForLife:
			self.offer = applyNoise(self.offerGene)
			#self.offer = random()			# xx JUST TESTING WITHOUT OFFER EVOLUTION !!
		else:
			self.offer = self.offerGene
		#Compute request
		if IsReactionNorm:
			self.request = self.reactionNorm.applyNorm(self.offer)


	def payoffFunction(self, myOffer, partnerOffer):
		self.payoff = partnerOffer - pow(myOffer, 2)
		#self.payoff = partnerOffer - 0.5 * pow(myOffer, 2)			# xx JUST TESTING
		#self.payoff = exp(partnerOffer * myOffer)	# xx JUST TESTING OTHER LOG-SUPERMODULAR PAYOFF FUNCTION


	def mutate(self):
		self.offerGene = mutateGene(self.offerGene)
		if self.offerGene < 0:
			self.offerGene = 0
		if IsReactionNorm:
			self.reactionNorm.mutate()
		else:
			self.request = mutateGene(self.request)



	@property
	def offerGene(self):
		return self._offerGene

	@offerGene.setter
	def offerGene(self, v):
       		self._offerGene = v

	@property
	def offer(self):
		return self._offer

	@offer.setter
	def offer(self, v):
       		self._offer = v

	@property
	def request(self):
		return self._request

	@request.setter
	def request(self, v):
       		self._request = v

	@property
	def partnerIndex(self):
		return self._partnerIndex

	@partnerIndex.setter
	def partnerIndex(self, v):
       		self._partnerIndex = v

	@property
	def payoff(self):
		return self._payoff

	@payoff.setter
	def payoff(self, v):
       		self._payoff = v

	@property
	def cumulatedPayoff(self):
		return self._cumulatedPayoff

	@cumulatedPayoff.setter
	def cumulatedPayoff(self, v):
       		self._cumulatedPayoff = v

	@property
	def fecundity(self):
		return self._fecundity

	@fecundity.setter
	def fecundity(self, v):
       		self._fecundity = v
