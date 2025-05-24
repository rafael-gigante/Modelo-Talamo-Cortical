# Noise factors
s2 = 1.5      # Additive white Gaussian noise strength
e2 = 0.5      # Threshold white Gaussian noise strength

# noise added to the membrane potential of the neurons
kisiSE = s2 * randn(nEs, nSim)
kisiME = s2 * randn(nEm, nSim)
kisiDE = s2 * randn(nEd, nSim)
kisiSI = s2 * randn(nINs, nSim)
kisiIret = s2 * randn(nIret, nSim)
kisiErel = s2 * randn(nErel, nSim)

# noise added to the threshold potential of the neurons
zetaSE = e2 * randn(nEs, nSim)
zetaME = e2 * randn(nEm, nSim)
zetaDE = e2 * randn(nEd, nSim)
zetaSI = e2 * randn(nINs, nSim)
zetaIret = e2 * randn(nIret, nSim)
zetaErel = e2 * randn(nErel, nSim)