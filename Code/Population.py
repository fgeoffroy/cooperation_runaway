from __future__ import division
import matplotlib.pyplot as plt
import sys
from random import *
from copy import *
from Individual import *
from Parameters import *
from Utilities import *
import numpy as np
import os



class Population:

	def __init__(self):
		self.individuals = []
		for i in range(PopulationSize):
			self.individuals.append(Individual())
		self.tempIndividuals = []
		self.solitaryIndexList = list(range(PopulationSize))
		self.pairedIndexList = []
		self._nbOfSolitaryIndividuals = PopulationSize
		self._nbOfPairedIndividuals = 0
		self._nbOfEncounters = 0
		self._nbOfInteractions = 0
		self._rateOfAcceptation = 0


	def evolve(self):
		if Plot:
			plt.ion()
			plt.figure("Generation")
			ax = plt.gca()
			ax.set_autoscale_on(True)
			hl, = ax.plot([], [],"-r")
			hl2, = ax.plot([], [],"-b")

		if File:
			if not os.path.exists("../Results"):
				os.mkdir("../Results")
			fileName = os.getcwd() + "/../Results/results" + repr(random()) + ".txt"
			with open(fileName, "wb") as results:
				results.write('"PopulationSize" "Lifespan" "NumberOfGenerations" "Beta" "Tau" "ProbabilityUniform" "MutationRate" "MutationStd" "NoiseForLife" "NoiseStd" "IsReactionNorm"')

				if IsReactionNorm:
					if IsDiscreteNorm:
						results.write(' "IsDiscreteNorm" "NbOfDiscreteRequests"')
					elif IsPolynomialNorm:
						results.write(' "IsPolynomialNorm" "NbOfPolynomialCoefficients"')
				results.write('\n')

				results.write(repr(PopulationSize) + " " + repr(Lifespan) + " " + repr(NumberOfGenerations) + " " + \
repr(Beta) + " " + repr(Tau) + " " + repr(ProbabilityUniform) + " " + repr(MutationRate) + " " + repr(MutationStd) + " " + \
repr(NoiseForLife) + " " + repr(NoiseStd) + " " + repr(IsReactionNorm))

				if IsReactionNorm:
					if IsDiscreteNorm:
						results.write(" " + repr(IsDiscreteNorm) + " " + repr(NbOfDiscreteRequests))
					elif IsPolynomialNorm:
						results.write(" " + repr(IsPolynomialNorm) + " " + repr(NbOfPolynomialCoefficients))
				results.write('\n')

				results.write('"generation" "offer" "request" "rateOfAcceptation"')
				if IsReactionNorm:
					if IsDiscreteNorm:
						for i in range(0, NbOfDiscreteRequests):
							results.write(' "request' + repr(i / NbOfDiscreteRequests) + '"')
					elif IsPolynomialNorm:
						for i in range(0, NbOfPolynomialCoefficients):
							results.write(' "a' + repr(i) + '"')
				results.write('\n')


		generation = 0
		while generation < NumberOfGenerations:
			self.resetGeneration()
			time = 0
			while time < Lifespan:
				threshold, timeUntilNextEvent = self.nextEvent()
				if (random() <= threshold):
					self.encounter()
				else:
					self.split()
				time += timeUntilNextEvent
			if self.nbOfEncounters != 0:
				self.rateOfAcceptation = self.nbOfInteractions / self.nbOfEncounters

			if generation % PrintStep == 0:
				if Print:
					print("generation :", generation)
					print("offer :", self.computeMeanOfferGene())
					print("request :", self.computeMeanRequest())
					print("")

				if Plot:
					update_line(hl, generation, self.computeMeanOfferGene())
					if IsReactionNorm:
						update_norm(self.computeMeanNorm())
					else:
						update_line(hl2, generation, self.computeMeanRequest())

				if File:
					with open(fileName, "a") as results:
						results.write(repr(generation) + " " + repr(self.computeMeanOfferGene()) + " " + repr(self.computeMeanRequest()) + " " + repr(self.rateOfAcceptation))
						if IsReactionNorm:
							if IsDiscreteNorm:
								meanRequestVector = self.computeMeanRequestVector()
								for req in meanRequestVector:
									results.write(" " + repr(req))
							elif IsPolynomialNorm:
								meanPolynomialCoeff = self.computeMeanPolynomialCoeff()
								for coeff in meanPolynomialCoeff:
									results.write(" " + repr(coeff))
						results.write("\n")

			self.reproduceAndDie()
			self.mutate()
			generation += 1


	def resetGeneration(self):
		self.nbOfEncounters = 0
		self.nbOfInteractions = 0
		self.solitaryIndexList = list(range(PopulationSize))
		self.pairedIndexList = []
		self.nbOfSolitaryIndividuals = PopulationSize
		self.nbOfPairedIndividuals = 0
		for i in self.individuals:
			i.reset()

	def nextEvent(self):
		lambdaExpo = (self.nbOfPairedIndividuals * Tau / 2) + self.nbOfSolitaryIndividuals * Beta
		timeUntilNextEvent = expovariate(lambdaExpo)
		for i in self.individuals:
			i.cumulatedPayoff = i.cumulatedPayoff + i.payoff * timeUntilNextEvent
		threshold = self.nbOfSolitaryIndividuals * Beta / lambdaExpo
		return (threshold, timeUntilNextEvent)


	def encounter(self):
		self.nbOfEncounters += 1
		#Pick the first individual
		i1 = randint(0, self.nbOfSolitaryIndividuals - 1)
		index1 = self.solitaryIndexList[i1]
		indiv1 = self.individuals[index1]
		self.solitaryIndexList.remove(index1)

		#Pick the second individual
		i2 = randint(0, self.nbOfSolitaryIndividuals - 2)
		index2 = self.solitaryIndexList[i2]
		indiv2 = self.individuals[index2]

		#Applying noise on offers if noise for interaction
		if NoiseForLife:
			offer1 = indiv1.offer
			offer2 = indiv2.offer
			request1 = indiv1.request
			request2 = indiv2.request
		else:
			offer1 = applyNoise(indiv1.offerGene)
			offer2 = applyNoise(indiv2.offerGene)
			if IsReactionNorm:
				request1 = indiv1.reactionNorm.applyNorm(offer1)
				request2 = indiv2.reactionNorm.applyNorm(offer2)
			else:
				request1 = indiv1.request
				request2 = indiv2.request

		if offer1 < request2 or offer2 < request1:
			self.solitaryIndexList.append(index1)
		else:
			self.nbOfInteractions += 1
			#Pairing individuals
			self.solitaryIndexList.remove(index2)
			self.pairedIndexList.append(index1)
			self.pairedIndexList.append(index2)
			self.nbOfSolitaryIndividuals = self.nbOfSolitaryIndividuals - 2
			self.nbOfPairedIndividuals = self.nbOfPairedIndividuals + 2
			indiv1.partnerIndex = index2
			indiv2.partnerIndex = index1
			#Collective action taking place
			indiv1.payoffFunction(offer1, offer2)
			indiv2.payoffFunction(offer2, offer1)


	def split(self):
		#Pick the first individual
		i1 = randint(0, self.nbOfPairedIndividuals - 1)
		index1 = self.pairedIndexList[i1]

		#Find his partner
		for i, elt in enumerate(self.pairedIndexList):
			if self.individuals[elt].partnerIndex == index1:
				i2 = i
				break
		index2 = self.pairedIndexList[i2]

		self.pairedIndexList.remove(index1)
		self.pairedIndexList.remove(index2)
		self.solitaryIndexList.append(index1)
		self.solitaryIndexList.append(index2)
		self.nbOfPairedIndividuals = self.nbOfPairedIndividuals - 2
		self.nbOfSolitaryIndividuals = self.nbOfSolitaryIndividuals + 2

		self.individuals[index1].partnerIndex = -1
		self.individuals[index2].partnerIndex = -1
		self.individuals[index1].payoff = 0
		self.individuals[index2].payoff = 0


	def computeMeanOfferGene(self):
		m = 0
		for i in self.individuals:
			m += i.offerGene
		m /= PopulationSize
		return m


	def computeMeanRequest(self):
		m = 0
		for i in self.individuals:
			m += i.request
		m /= PopulationSize
		return m


	def computeMeanNorm(self):
		m = np.empty(shape=(0,0))
		for x in np.linspace(0, 1, PlotNormNbOfSteps):
			y = 0
			for i in self.individuals:
				y += i.reactionNorm.applyNorm(x)
			y /= PopulationSize
			m = np.append(m, y)
		return m


	def computeMeanRequestVector(self):
		m = np.zeros(NbOfDiscreteRequests)
		for i in self.individuals:
			m += i.reactionNorm.requestVector
		m = m / PopulationSize
		return m


	def computeMeanPolynomialCoeff(self):
		m = np.zeros(NbOfPolynomialCoefficients)
		for i in self.individuals:
			m += i.reactionNorm.polynomeCoefficients
		m = m / PopulationSize
		return m


	def reproduceAndDie(self):
		#Normalization (min fecundity = 0)
		minFecundity = sys.maxint
		for i in self.individuals:
			if i.cumulatedPayoff < minFecundity:
				minFecundity = i.cumulatedPayoff
		for i in self.individuals:
			i.fecundity = i.cumulatedPayoff - minFecundity
		#Normalization (relative fecundity)
		totalFecundity = 0
		for i in self.individuals:
			totalFecundity += i.fecundity
		if totalFecundity != 0:
			for i in self.individuals:
				i.fecundity = i.fecundity / totalFecundity
		#Switch to tempIndividuals
		self.tempIndividuals = deepcopy(self.individuals)
		#Reproduce stochastically accordingly to relative fecundity (or randomly if same fecundity)
		if totalFecundity == 0:
			i = 0
			while i < PopulationSize:
				self.individuals[i] = deepcopy(choice(self.tempIndividuals))
				i += 1
		else:
			i = 0
			while i < PopulationSize:
				r = random()
				comp = 0
				shuffle(self.tempIndividuals)
				j = 0
				while j < PopulationSize:
					comp += self.tempIndividuals[j].fecundity
					if comp >= r:
						self.individuals[i] = deepcopy(self.tempIndividuals[j])
						break
					j += 1
				i += 1


	def mutate(self):
		for i in self.individuals:
			i.mutate()


	@property
	def nbOfSolitaryIndividuals(self):
		return self._nbOfSolitaryIndividuals

	@nbOfSolitaryIndividuals.setter
	def nbOfSolitaryIndividuals(self, v):
       		self._nbOfSolitaryIndividuals = v

	@property
	def nbOfPairedIndividuals(self):
		return self._nbOfPairedIndividuals

	@nbOfPairedIndividuals.setter
	def nbOfPairedIndividuals(self, v):
       		self._nbOfPairedIndividuals = v

	@property
	def nbOfEncounters(self):
		return self._nbOfEncounters

	@nbOfEncounters.setter
	def nbOfEncounters(self, v):
       		self._nbOfEncounters = v

	@property
	def nbOfInteractions(self):
		return self._nbOfInteractions

	@nbOfInteractions.setter
	def nbOfInteractions(self, v):
       		self._nbOfInteractions = v

	@property
	def rateOfAcceptation(self):
		return self._rateOfAcceptation

	@rateOfAcceptation.setter
	def rateOfAcceptation(self, v):
       		self._rateOfAcceptation = v
