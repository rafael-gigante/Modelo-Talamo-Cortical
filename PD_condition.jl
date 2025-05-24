# Division Factor:
fac = fac_PD

# Factor for coupling strengths
r_s = rand(nEs) # Layer S
r_m = rand(nEm) # Layer M
r_d = rand(nEd) # Layer D
r_ins = rand(nINs) # INs
r_ret = rand(nIret) # Reticular nucleus
r_rel = rand(nErel) # Relay nucleus

# COUPLING STRENGTHS within each structure
W_EEs = (-5e1 / fac) .* r_s # S to S
W_EEm = (-5e1 / fac) .* r_m # M to M
W_EEd = (-5e1 / fac) .* r_d # D to D
W_IIins = (-5e1 / fac) .* r_ins # INs to INs
W_IIret = (-5e1 / fac) .* r_ret # Reticular to Reticular
W_EErel = (0 / fac) .* r_rel # Relay to Relay

# COUPLING STRENGTHS between structures (PD)
# S
W_EEsm = (3e2 / fac) .* r_s # M to S
W_EEsd = (5e2 / fac) .* r_s # D to S
W_EIsINs = (-7.5e2 / fac) .* r_s # INs to S
W_EIsRet = (0 / fac) .* r_s # Reticular to S
W_EEsRel = (0 / fac) .* r_s # Relay to S

# M
W_EEms = (1e1 / fac) .* r_m # S to M
W_EEmd = (0 / fac) .* r_m # D to M
W_EImINs = (-7.5e2 / fac) .* r_m # INs to M
W_EImRet = (0 / fac) .* r_m # Reticular to M
W_EEmRel = (0 / fac) .* r_m # Relay to M

# D
W_EEds = (3e2 / fac) .* r_d # S to D
W_EEdm = (0 / fac) .* r_d # M to D
W_EIdINs = (-5e3 / fac) .* r_d # INs to D
W_EIdRet = (0 / fac) .* r_d # Reticular to D
W_EEdRel = (1e3 / fac) .* r_d # Relay to D

# INs
W_IE_INs_s = (2e2 / fac) .* r_ins # S to INs
W_IE_INs_m = (2e2 / fac) .* r_ins # M to INs
W_IE_INs_d = (2e2 / fac) .* r_ins # D to INs
W_II_INs_Ret = (0 / fac) .* r_ins # Reticular to INs
W_IE_INs_Rel = (1e3 / fac) .* r_ins # Relay to INs

# Ret.
W_IE_Ret_s = (0 / fac) .* r_ret # S to Reticular
W_IE_Ret_m = (0 / fac) .* r_ret # M to Reticular
W_IE_Ret_d = (1e2 / fac) .* r_ret # D to Reticular
W_II_Ret_INs = (0 / fac) .* r_ret # INs to Reticular
W_IE_Ret_Rel = (5e2 / fac) .* r_ret # Relay to Reticular

# Rel.
W_EERels = (0 / fac) .* r_rel # S to Relay
W_EERelm = (0 / fac) .* r_rel # M to Relay
W_EEReld = (1e2 / fac) .* r_rel # D to Relay
W_EIRelINs = (0 / fac) .* r_rel # INs to Relay
W_EIRelRet = (-2.5e3 / fac) .* r_rel # Reticular to Relay

# Construct the normalized mean synaptic weights
# 1=Layer S, 2=Layer M, 3=Layer D, 4=INs, 5=TCR, 6=TRN
Z = zeros(6,6)  # Initialize a 6x6 matrix with zeros

Z[1,1] = mean(W_EEs)
Z[2,2] = mean(W_EEm)
Z[3,3] = mean(W_EEd)
Z[4,4] = mean(W_IIins)
Z[6,6] = mean(W_IIret)
Z[5,5] = mean(W_EErel)

Z[2,1] = mean(W_EEsm)
Z[3,1] = mean(W_EEsd)
Z[4,1] = mean(W_EIsINs)
Z[6,1] = mean(W_EIsRet)
Z[5,1] = mean(W_EEsRel)

Z[1,2] = mean(W_EEms)
Z[3,2] = mean(W_EEmd)
Z[4,2] = mean(W_EImINs)
Z[6,2] = mean(W_EImRet)
Z[5,2] = mean(W_EEmRel)

Z[1,3] = mean(W_EEds)
Z[2,3] = mean(W_EEdm)
Z[4,3] = mean(W_EIdINs)
Z[6,3] = mean(W_EIdRet)
Z[5,3] = mean(W_EEdRel)

Z[1,4] = mean(W_IE_INs_s)
Z[2,4] = mean(W_IE_INs_m)
Z[3,4] = mean(W_IE_INs_d)
Z[6,4] = mean(W_II_INs_Ret)
Z[5,4] = mean(W_IE_INs_Rel)

Z[1,6] = mean(W_IE_Ret_s)
Z[2,6] = mean(W_IE_Ret_m)
Z[3,6] = mean(W_IE_Ret_d)
Z[4,6] = mean(W_II_Ret_INs)
Z[5,6] = mean(W_IE_Ret_Rel)

Z[1,5] = mean(W_EERels)
Z[2,5] = mean(W_EERelm)
Z[3,5] = mean(W_EEReld)
Z[4,5] = mean(W_EIRelINs)
Z[6,5] = mean(W_EIRelRet)

ZZ_pd = copy(Z)
Zx_pd = Z / maximum(Z)  # Normalize by the max value in Z