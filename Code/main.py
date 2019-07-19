from __future__ import division
from Population import *
import sys

if IsReactionNorm:
	count = 0
	if IsDiscreteNorm:
		count += 1
	if IsPolynomialNorm:
		count += 1
	if count > 1:
		sys.exit("Error: Polynomial and discrete norms are allowed at the same time")
	if count == 0:
		sys.exit("Error: No method for the reaction norm is indicated")

population = Population()
population.evolve()
