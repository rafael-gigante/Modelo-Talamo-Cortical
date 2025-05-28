using Statistics
using DelimitedFiles
using Random
###########################################################################
# ECP ligada ou desligada:
fidD = 5 * 67 # amplitude da corrente de DBS
#fidD = 0

sim_time = 15   # tempo total de simulação em segundos (deve ser múltiplo de 3 se ECP ligada)
folder = "dbs_n"
trial = 10
###########################################################################
include("simulation_params.jl")
###########################################################################
# Neurônios:
# Define o número de neurônios e a distribuição de tipos deles em cada estrutura e seus parâmetros
# com variações neuronais seguindo o algoritmo de Izhikevich
include("neurons.jl")
###########################################################################
# SYNAPSES COUPLING MATRIX:
fac_N  = 2.5 
fac_PD = 5

include("Normal_condition.jl")
include("PD_condition.jl")

ZZ_diff = ZZ_pd - ZZ_normal
Zx_diff =  ZZ_diff / (maximum(ZZ_diff)) # Normalização

# Salvando a matriz de acoplamento da diferença entre as condições normal e PD
writedlm("$(folder)/trial_$(trial)/ZZ_diff.csv", ZZ_diff, ',')
writedlm("$(folder)/trial_$(trial)/Zx_diff.csv", Zx_diff, ',')
###########################################################################
# Ruido branco aditivo e limiar:
include("white_noise.jl")
###########################################################################
# Atividade de fundo poissoniana:
include("poissonspikes.jl")
###########################################################################
# Parametros da ECP:
include("dbs.jl")
writedlm("$(folder)/trial_$(trial)/I_dbs_pre.csv", I_dbs[1, :], ',')
writedlm("$(folder)/trial_$(trial)/I_dbs_post.csv", I_dbs[2, :], ',')
###########################################################################
vr = -65.0  # Valor de repouso da membrana

# Potenciais de membrana iniciais e variáveis de recuperação:
vEs, uEs = fill(vr, nEs, nSim), zeros(nEs, nSim)
vEm, uEm = fill(vr, nEm, nSim), zeros(nEm, nSim)
vEd, uEd = fill(vr, nEd, nSim), zeros(nEd, nSim)
vIs, uIs = fill(vr, nINs, nSim), zeros(nINs, nSim)
vIret, uIret = fill(vr, nIret, nSim), zeros(nIret, nSim)
vErel, uRel = fill(vr, nErel, nSim), zeros(nErel, nSim)

# Correntes sinápticas
EPSCs = zeros(nSim)
EPSCm = zeros(nSim)
EPSCd = zeros(nSim)
IPSC_INs = zeros(nSim)
EPSC_rel = zeros(nSim)
IPSC_ret = zeros(nSim)

EPSCdF = zeros(nSim)
EPSC_relD = zeros(nSim)

# Valor inicial das váriaveis de recuperação e sinápticas:
rEs, xEs, IsEs = zeros(3), ones(3), zeros(3)
rEm, xEm, IsEm = zeros(3), ones(3), zeros(3)
rEd, xEd, IsEd = zeros(3), ones(3), zeros(3)
rINs, xINs, IsINs = zeros(3), ones(3), zeros(3)
rIret, xIret, IsIret = zeros(3), ones(3), zeros(3)
rErel, xErel, IsErel = zeros(3), ones(3), zeros(3)

rEdF, xEdF, IsEdF = zeros(3), ones(3), zeros(3)
rErelD, xErelD, IsErelD = zeros(3), ones(3), zeros(3)
###########################################################################
# Inclui as estruturas corticais e tálamicas:
include("cortex.jl")
include("thalamus.jl")
using .Cortex
using .Thalamus

for i in tVec
    # Camada cortical S
    global rEs, xEs, IsEs
    vEs[:, i+1], uEs[:, i+1], rEs, xEs, IsEs, EPSCs[i+1] =
    S_layer(aEs, bEs, cEs, dEs, nEs, vEs[:, i], uEs[:, i], rEs, xEs, IsEs,
    EPSCs[i-τ_wL-τ_syn], EPSCd[i-τ_L-τ_syn], EPSCm[i-τ_L-τ_syn], EPSC_rel[i-τ_TC-τ_syn], IPSC_INs[i-τ_wL-τ_syn], IPSC_ret[i-τ_TC-τ_syn],
    W_EEs, W_EEsd, W_EEsm, W_EEsRel, W_EIsINs, W_EIsRet,
    I_ps[1,1,i-τ_wL-τ_syn], I_ps[1,2,i-τ_wL-τ_syn], kisiSE[:,i], zetaSE[:,i], IdcS_E, fidS .* I_dbs[2,i], n_conn_S, dt)
    
    # Camada cortical M
    global rEm, xEm, IsEm
    vEm[:, i+1], uEm[:, i+1], rEm, xEm, IsEm, EPSCm[i+1] =
    M_layer(aEm, bEm, cEm, dEm, nEm, vEm[:, i], uEm[:, i], rEm, xEm, IsEm,
    EPSCm[i-τ_wL-τ_syn], EPSCs[i-τ_L-τ_syn], EPSCd[i-τ_L-τ_syn], EPSC_rel[i-τ_TC-τ_syn], IPSC_INs[i-τ_wL-τ_syn], IPSC_ret[i-τ_TC-τ_syn],
    W_EEm, W_EEms, W_EEmd, W_EEmRel, W_EImINs, W_EImRet,
    I_ps[2,1,i-τ_wL-τ_syn], I_ps[2,2,i-τ_wL-τ_syn], kisiME[:,i], zetaME[:,i], IdcM_E, fidM .* I_dbs[2,i], n_conn_M, dt)
    
    # Camada cortical P
    global rEd, xEd, IsEd, rEdF, xEdF, IsEdF
    vEd[:, i+1], uEd[:, i+1], rEd, xEd, IsEd, EPSCd[i+1], rEdF, xEdF, IsEd, EPSCdF[i+1] =
    D_layer(aEd, bEd, cEd, dEd, nEd, n_hyp, vEd[:, i], uEd[:, i], rEd, xEd, IsEd, rEdF, xEdF, IsEdF,
    EPSCd[i-τ_wL-τ_syn], EPSCs[i-τ_L-τ_syn], EPSCm[i-τ_L-τ_syn], EPSC_relD[i-τ_TC-τ_syn], IPSC_INs[i-τ_wL-τ_syn], IPSC_ret[i-τ_TC-τ_syn],
    W_EEd, W_EEds, W_EEdm, W_EEdRel, W_EIdINs, W_EIdRet,
    I_ps[3,1,i-τ_wL-τ_syn], I_ps[3,2,i-τ_wL-τ_syn], kisiDE[:,i], zetaDE[:,i], IdcD_E, fidD .* I_dbs[:,i], dt)
    
    # Interneurônios corticais
    global rINs, xINs, IsINs
    vIs[:, i+1], uIs[:, i+1], rINs, xINs, IsINs, IPSC_INs[i+1] =
    ctx_INs(aIs, bIs, cIs, dIs, nINs, vIs[:, i], uIs[:, i], rINs, xINs, IsINs,
    IPSC_INs[i-τ_wL-τ_syn], EPSCs[i-τ_L-τ_syn], EPSCm[i-τ_L-τ_syn], EPSCd[i-τ_L-τ_syn], EPSC_rel[i-τ_TC-τ_syn], IPSC_ret[i-τ_TC-τ_syn],
    W_IE_INs_d, W_IE_INs_s, W_IE_INs_m, W_II_INs_Ret, W_IIins, W_IE_INs_Rel,
    I_ps[4,1,i-τ_wL-τ_syn], I_ps[4,2,i-τ_wL-τ_syn], kisiSI[:,i], zetaSI[:,i], Idc_INs, fidCI .* I_dbs[2,i], n_conn_CI, dt)
    
    # Neurônios reticulares talâmicos
    global rIret, xIret, IsIret
    vIret[:, i+1], uIret[:, i+1], rIret, xIret, IsIret, IPSC_ret[i+1] =
    thm_ret(aIret, bIret, cIret, dIret, nIret, vIret[:, i], uIret[:, i], rIret, xIret, IsIret,
    IPSC_ret[i-τ_wL-τ_syn], EPSCs[i-τ_CT-τ_syn], EPSCm[i-τ_CT-τ_syn], EPSCdF[i-τ_CT-τ_syn], IPSC_INs[i-τ_CT-τ_syn], EPSC_rel[i-τ_L-τ_syn],
    W_IIret, W_IE_Ret_s, W_IE_Ret_m, W_IE_Ret_d, W_II_Ret_INs, W_IE_Ret_Rel,
    0 .* I_ps[5,1,i-τ_wL-τ_syn], 0 .* I_ps[5,2,i-τ_wL-τ_syn], kisiIret[:,i], zetaIret[:,i], Idc_Ret, fidN .* I_dbs[2,i], n_conn_N, dt)
    
    # Neurônios de relé talâmicos
    global rErel, xErel, IsErel, rErelD, xErelD, IsErelD
    vErel[:, i+1], uRel[:, i+1], rErel, xErel, IsErel, EPSC_rel[i+1], rErelD, xErelD, IsErelD, EPSC_relD[i+1] =
    thm_rel(aErel, bErel, cErel, dErel, nErel, vErel[:, i], uRel[:, i], rErel, xErel, IsErel, rErelD, xErelD, IsErelD,
    EPSC_rel[i-τ_wL-τ_syn], EPSCs[i-τ_CT-τ_syn], EPSCm[i-τ_CT-τ_syn], EPSCdF[i-τ_CT-τ_syn], IPSC_INs[i-τ_CT-τ_syn], IPSC_ret[i-τ_L-τ_syn],
    W_EErel, W_EERels, W_EERelm, W_EEReld, W_EIRelINs, W_EIRelRet,
    0 .* I_ps[6,1,i-τ_wL-τ_syn], 0 .* I_ps[6,2,i-τ_wL-τ_syn], kisiErel[:,i], zetaErel[:,i], Idc_Rel, fidR .* I_dbs[2,i], n_conn_R, dt)
    println("Simulation step: ", i, " / ", nSim)
    
end

vEs = vEs[:, chop_till+1:nSim]
vEm = vEm[:, chop_till+1:nSim]
vEd = vEd[:, chop_till+1:nSim]
vIs = vIs[:, chop_till+1:nSim]
vIret = vIret[:, chop_till+1:nSim]
vErel = vErel[:, chop_till+1:nSim]

EPSCs = EPSCs[chop_till+1:nSim]
EPSCm = EPSCm[chop_till+1:nSim]
EPSCd = EPSCd[chop_till+1:nSim]
IPSC_INs = IPSC_INs[chop_till+1:nSim]
IPSC_ret = IPSC_ret[chop_till+1:nSim]
EPSC_rel = EPSC_rel[chop_till+1:nSim]
EPSCdF = EPSCdF[chop_till+1:nSim]
EPSC_relD = EPSC_relD[chop_till+1:nSim]

# Cálculo do LFP
ρ = 0.27
r = 100e-6
LFP = (EPSCd - IPSC_INs)/(4 * pi * ρ * r)

writedlm("$(folder)/trial_$(trial)/vEs.csv", vEs, ',')
writedlm("$(folder)/trial_$(trial)/vEm.csv", vEm, ',')
writedlm("$(folder)/trial_$(trial)/vEd.csv", vEd, ',')
writedlm("$(folder)/trial_$(trial)/vIs.csv", vIs, ',')
writedlm("$(folder)/trial_$(trial)/vIret.csv", vIret, ',')
writedlm("$(folder)/trial_$(trial)/vErel.csv", vErel, ',')
writedlm("$(folder)/trial_$(trial)/cEs.csv", cEs, ',')
writedlm("$(folder)/trial_$(trial)/cEm.csv", cEm, ',')
writedlm("$(folder)/trial_$(trial)/cEd.csv", cEd, ',')
writedlm("$(folder)/trial_$(trial)/cIs.csv", cIs, ',')
writedlm("$(folder)/trial_$(trial)/cIret.csv", cIret, ',')
writedlm("$(folder)/trial_$(trial)/cErel.csv", cErel, ',')
writedlm("$(folder)/trial_$(trial)/LFP.csv", LFP, ',')