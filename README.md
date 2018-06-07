# BayesianAxisAlliesCalculator


This program calculates combat outcome probabilities for Axis and Allies battles. Most of the available calculators online do Monte Carlo simulation, where each battle is simulated (using a random number generator in place of dice rolls) many times, and a count of each result is listed. In contrast, this work takes a bayesian approach, calculating the absolute probability of a node given the probability of its descendence from its parents and the absolute probability of its parents.
