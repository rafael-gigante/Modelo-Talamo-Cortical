# Fatores de escala para a força do ruído aditivo branco e do limiar do ruído branco
s2 = 1.5      # Força do ruído aditivo branco
e2 = 0.5      # Força do limiar do ruído branco

# Ruído adicionado ao potencial de membrana dos neurônios
kisiSE = s2 * randn(nEs, nSim)
kisiME = s2 * randn(nEm, nSim)
kisiDE = s2 * randn(nEd, nSim)
kisiSI = s2 * randn(nINs, nSim)
kisiIret = s2 * randn(nIret, nSim)
kisiErel = s2 * randn(nErel, nSim)

# Ruído adicionado ao potencial de limiar dos neurônios
zetaSE = e2 * randn(nEs, nSim)
zetaME = e2 * randn(nEm, nSim)
zetaDE = e2 * randn(nEd, nSim)
zetaSI = e2 * randn(nINs, nSim)
zetaIret = e2 * randn(nIret, nSim)
zetaErel = e2 * randn(nErel, nSim)