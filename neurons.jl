function neuron_parameters(a, b, c, d, n1, n2, T1, T2, neuron_type)
    if neuron_type == "E"
        re1, re2 = rand(n1), rand(n2)
        aE = vcat(fill(a[T1], n1), fill(a[T2], n2))
        bE = vcat(fill(b[T1], n1), fill(b[T2], n2))
        cE = vcat(c[T1] .+ (15 .* re1 .^ 2), c[T2] .+ (15 .* re2 .^ 2))
        dE = vcat(d[T1] .- (0.6 .* re1 .^ 2), d[T2] .- (0.6 .* re2 .^ 2))
        return aE, bE, cE, dE
    elseif neuron_type == "I"
        ri1, ri2 = rand(n1), rand(n2)
        aI = vcat(a[T1] .+ (0.008 .* ri1), a[T2] .+ (0.008 .* ri2))
        bI = vcat(b[T1] .- (0.005 .* ri1), b[T2] .- (0.005 .* ri2))
        cI = vcat(fill(c[T1], n1), fill(c[T2], n2))
        dI = vcat(fill(d[T1], n1), fill(d[T2], n2))
        return aI, bI, cI, dI
    end
end

# Número de neurônios em cada estrutura
# S layer, M layer, D layer, INs, TCR, TRN
nEs, nEm, nEd, nErel, nINs, nIret = 100, 100, 100, 100, 100, 40
n_tot = nEs + nEm + nEm + nINs + nErel + nIret

# Distribuição de neurônios em cada estrutura
nE1s, nE2s = Int(0.5 * nEs), Int(0.5 * nEs)
nE1m, nE2m = Int(1 * nEm), Int(0 * nEm)
nE1d, nE2d = Int(0.7 * nEd), Int(0.3 * nEd)
nErel1, nErel2 = Int(0.7 * nErel), Int(0.3 * nErel)
nINs1, nINs2 = Int(0.5 * nINs), Int(0.5 * nINs)
nIret1, nIret2 = Int(0.5 * nIret), Int(0.5 * nIret)

# Tipos de neurônios excitatórios e inibitórios em cada estrutura
E1s, E2s = 1, 2
E1m, E2m = 1, 1
E1d, E2d = 1, 2
I1s, I2s = 3, 4
I1ret, I2ret = 6, 6
E1rel, E2rel = 5, 5

# Parâmetros dos neurônios
# 1)RS  2)IB   3)FS   4)LTS   5)Rel(TC)   6)Ret(TR)
a = [0.02, 0.02, 0.1, 0.02, 0.02, 0.02]
b = [0.2, 0.2, 0.2, 0.25, 0.25, 0.25]
c = [-65, -55, -65, -65, -65, -65]
d = [8, 4, 2, 2, 0.05, 2.05]

# Diversidade dos parâmetros dos neurônios de acordo com o algoritmo de Izhikevich
# Camada S
aEs, bEs, cEs, dEs = neuron_parameters(a, b, c, d, nE1s, nE2s, E1s, E2s, "E")
# Camada M
aEm, bEm, cEm, dEm = neuron_parameters(a, b, c, d, nE1m, nE2m, E1m, E2m, "E")
# Camada P
aEd, bEd, cEd, dEd = neuron_parameters(a, b, c, d, nE1d, nE2d, E1d, E2d, "E")
# Neurônios Relé
aErel, bErel, cErel, dErel = neuron_parameters(a, b, c, d, nErel1, nErel2, E1rel, E2rel, "I")
# INs
aIs, bIs, cIs, dIs = neuron_parameters(a, b, c, d, nINs1, nINs2, I1s, I2s, "I")
# Neurônios Reticulares
aIret, bIret, cIret, dIret = neuron_parameters(a, b, c, d, nIret1, nIret2, I1ret, I2ret, "I")
###########################################################################
# Corrente DC:
# 1)RS  2)IB   3)FS   4)LTS   5)Rel(TC)   6)Ret(TR)
Idc = [3.5, 3.6, 3.8, 0.4, 0.6, 0.6] .+ 0.1

IdcS_E = vcat(fill(Idc[E1s], nE1s), fill(Idc[E2s], nE2s))
IdcM_E = vcat(fill(Idc[E1m], nE1m), fill(Idc[E2m], nE2m))
IdcD_E = vcat(fill(Idc[E1d], nE1d), fill(Idc[E2d], nE2d))
Idc_INs = vcat(fill(Idc[I1s], nINs1), fill(Idc[I2s], nINs2))
Idc_Ret = vcat(fill(Idc[I1ret], nIret1), fill(Idc[I2ret], nIret2))
Idc_Rel = vcat(fill(Idc[E1rel], nErel1), fill(Idc[E2rel], nErel2))