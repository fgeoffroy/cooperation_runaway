from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from random import *
from copy import *
from Parameters import *


def update_line(hl, generation, value):
	plt.ion()
	plt.figure("Generation")
	hl.set_xdata(np.append(hl.get_xdata(), generation))
	hl.set_ydata(np.append(hl.get_ydata(), value))
	ax = plt.gca()
	ax.relim()
	ax.autoscale_view(True,True,True)
	plt.draw()
	plt.pause(0.00001)


def update_norm(requestList):
	plt.figure("Norm")
	plt.clf()
	plt.xlim([0,1])
	plt.ylim([0,1])
	plt.plot(np.linspace(0, 1, PlotNormNbOfSteps), requestList)
	plt.draw()
	plt.pause(0.00001)


def applyNoise(offerGene):
	x = gauss(offerGene, NoiseStd)
	if x >= 0:
		return x
	else:
		return 0


def mutateGene(gene):
	if random() <= MutationRate:
		if random() <= ProbabilityUniform:
			#return uniform(0, 1)
			return random()
		else:
			return gauss(gene, MutationStd)
	else:
		return gene


def mutateWeight(weight):
	if random() <= MutationRate:
		if random() <= ProbabilityUniform:
			return uniform(LowerWeightLimit, UpperWeightLimit)
		else:
			weight = gauss(weight, MutationStd)
			if weight > UpperWeightLimit:
				return UpperWeightLimit
			elif weight < LowerWeightLimit:
				return LowerWeightLimit
			else:
				return weight
	else:
		return weight


vMutateWeight = np.vectorize(mutateWeight)


vMutateDiscreteNorm = np.vectorize(mutateGene)


vMutatePolynomialNorm = np.vectorize(mutateGene)
