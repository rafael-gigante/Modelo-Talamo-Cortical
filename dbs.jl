function dbs_delta(fdbs, dbs_duration, dev, nSim, Fs, cut)
    Tdbs = 1 / fdbs  # intervalo entre pulsos de DBS (em segundos)
    N_steps = Fs * Tdbs  # número de passos entre os pulsos de DBS (Fs = 1000/dt Hz)
    dbs = 1:round(Int, N_steps):dbs_duration  # Índices dos pulsos de DBS

    I_dbs_full = zeros(dbs_duration)
    I_dbs_full[dbs] .= 1  # define os pulsos de DBS

    if dev == 1
        return I_dbs_full
    else
        return vcat(
            zeros(Int32(cut)),  # Preenchimento inicial com zeros
            zeros(Int32((nSim - cut) ÷ dev)),  # Primeiro segmento inativo
            I_dbs_full,  # Segmento ativo de DBS
            zeros(Int32((nSim - cut) ÷ dev))  # Último segmento inativo
        )
    end
end

function TMsynE_dbs(dbs, nSim, tdelay, dt=0.1)
    # parametros sinápticos para sinapses facilitadoras, depressoras e pseudo-lineares
    τf = [670, 17, 326]    # Constantes de tempo de facilitação (ms)
    τd = [138, 671, 329]   # Constantes de tempo de depressão (ms)
    U = [0.09, 0.5, 0.29]    # Fator de utilização (probabilidade de liberação)
    A = [0.20, 0.63, 0.17]   # Eficácia sináptica
    τsE = 3               # Constante de tempo de decaimento da corrente sináptica (ms)

    # Inicializa as variáveis de estado
    r = zeros(3, nSim)  # recuperação de neurotransmissores
    x = ones(3, nSim)   # recursos de neurotransmissores disponíveis
    Is = zeros(3, nSim) # Contribuições da corrente sináptica

    # Itera sobre os tipos de sinapse (F, D, P)
    for p in 1:3
        for i in (tdelay + 1):(nSim - 1)
            r[p, i + 1] = r[p, i] + dt * (-r[p, i] / τf[p] + U[p] * (1 - r[p, i]) * dbs[i - tdelay])
            x[p, i + 1] = x[p, i] + dt * ((1 / τd[p]) * (1 - x[p, i]) - r[p, i + 1] * x[p, i] * dbs[i - tdelay])
            Is[p, i + 1] = Is[p, i] + dt * ((-Is[p, i] / τsE) + A[p] * r[p, i + 1] * x[p, i] * dbs[i - tdelay])
        end
    end

    # Corrente pós-sináptica total devido a ECP
    I_dbs = sum(Is, dims=1)[1, :]

    return I_dbs
end

# Define the % of affected neurons in D layer by DBS
nh = 0.1
# Impacto da ECP nas estruturas neuronais conectadas à camada D
fidCI = abs(1 * fidD)   # Fidelidade sináptica para portadores de ECP invadirem CIs
fidM = abs(0 * fidD)    # Fidelidade sináptica para portadores de ECP invadirem a camada M
fidS = abs(1 * fidD)    # Fidelidade sináptica para portadores de ECP invadirem a camada S
fidR = abs(1 * fidD)    # Fidelidade sináptica para portadores de ECP invadirem a camada TCR
fidN = abs(1 * fidD)    # Fidelidade sináptica para portadores de ECP invadirem a camada TRN

nCI, nS, nR, nN = 1, 1, 1, 1
n_conn_CI = nCI * nh * nINs     # Percentual de neurônios CI com contato sináptico das arvores axonais hiperdiretas
n_conn_S = nS * nh * nEs        # Percentual de neurônios S com contato sináptico das arvores axonais hiperdiretas
n_conn_M = 0 * nh * nEm         # Percentual de neurônios M com contato sináptico das arvores axonais hiperdiretas
n_conn_R = nR * nh * nErel      # Percentual de neurônios R com contato sináptico das arvores axonais hiperdiretas
n_conn_N = nN * nh * nIret      # Percentual de neurônios N com contato sináptico das arvores axonais hiperdiretas

n_hyp = nEd * nh  # Número de neurônios hiperdiretos

# Inicializa a corrente de DBS
I_dbs = zeros(2, nSim)
if fidD != 0
    dev = 3
    fdbs = 130

    dbs_duration = (nSim - chop_till) ÷ dev
    I_dbs_pre = dbs_delta(fdbs, Int32(dbs_duration), dev, nSim, Fs, chop_till)

    # Pulsos de DBS pós-sinápticos (intra-axonal)
    I_dbs_post = TMsynE_dbs(I_dbs_pre, nSim, τ_syn, dt)

    I_dbs[1, :] .= I_dbs_pre
    I_dbs[2, :] .= I_dbs_post
end