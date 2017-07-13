import numpy
from mcmc import mcmc_framework
from mcmc_examples.lighthouse.lighthouse_model import LightHouseModel
from scipy.io import wavfile

# Data to use for 10 data points
data_points = [
    2.37706, 4.51142, 0.406605, 0.909418, 0.642899,
    1.21925, 1.47647, -2.95771, -0.801802, -1.86529
]


# Run and get the samples
model = LightHouseModel(
    alpha_jump_scale, alpha_min, alpha_max, beta_jump_scale, beta_min,
    beta_max)
samples = mcmc_framework.run_mcmc(
    model, data_points, n_samples, seed=seed, n_chips=23*48)

# Save the results
numpy.save("results.npy", samples)
numpy.savetxt("results.csv", samples, fmt="%f", delimiter=",")
