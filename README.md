continuous_time_ca_sampler
==========================

The code takes as an input a time series vector of calcium observations
and produces samples from the posterior distribution of the underlying
spike in continuous time. The code also samples the model parameters
(baseline, spike amplitude, initial calcium concentration, firing rate,
noise variance) and also iteratively re-estimates the discrete time
constant of the model. More info can be found at

Pnevmatikakis, E., Merel, J., Pakman, A. &amp; Paninski, L. (2014).
Bayesian spike inference from calcium imaging data. Asilomar Conf. on
Signals, Systems, and Computers. http://arxiv.org/abs/1311.6864

