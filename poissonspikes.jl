function poisson_spike_gen(fr, tSim, nTrials, dt)
    nBins = Int(floor(tSim / dt))  # número de bins
    spikeMat = rand(nTrials, nBins) .< (fr * dt)  # Processo de Poisson
    return spikeMat
end

function TMsynE(tevent, nSim, tdelay, dt)
    τf = [670, 17, 326]
    τd = [138, 671, 329]
    U = [0.09, 0.5, 0.29]
    A = [0.20, 0.63, 0.17]  

    τsE = 3  # Constante de tempo de decaimento da corrente sináptica

    # inicializa as variáveis de estado
    r = zeros(3, nSim)
    x = ones(3, nSim)
    Is = zeros(3, nSim)
    spd = zeros(nSim)

    # Certifique-se de que tevent esteja dentro dos limites
    tevent = tevent[tevent .<= nSim]
    spd[tevent] .= 1 / dt

    # Itera sobre os tipos de sinapse (Facilitadora, Depressora, Pseudo-linear)
    for p in 1:3
        for i in (tdelay + 1):(nSim - 1)
            r[p, i + 1] = r[p, i] + dt * (-r[p, i] / τf[p] + U[p] * (1 - r[p, i]) * spd[i - tdelay])
            x[p, i + 1] = x[p, i] + dt * ((1 / τd[p]) * (1 - x[p, i]) - r[p, i + 1] * x[p, i] * spd[i - tdelay])
            Is[p, i + 1] = Is[p, i] + dt * ((-Is[p, i] / τsE) + A[p] * r[p, i + 1] * x[p, i] * spd[i - tdelay])
        end
    end

    return vec(sum(Is, dims=1))
end

function TMsynI(tevent, nSim, tdelay, dt)
    τf = [376, 21, 62]
    τd = [45, 706, 144]
    U = [0.016, 0.25, 0.32]
    A = [0.08, 0.75, 0.17]  

    τsI = 11  # decaimento da corrente sináptica

    # Inicializa as variáveis de estado
    r = zeros(3, nSim)
    x = ones(3, nSim)
    Is = zeros(3, nSim)
    spd = zeros(nSim)

    # Certifique-se de que tevent esteja dentro dos limites 
    tevent = tevent[tevent .<= nSim]
    spd[tevent] .= 1 / dt  

    # Itera sobre os tipos de sinapse (Facilitadora, Depressora, Pseudo-linear)
    for p in 1:3
        for i in (tdelay + 1):(nSim - 1)
            r[p, i + 1] = r[p, i] + dt * (-r[p, i] / τf[p] + U[p] * (1 - r[p, i]) * spd[i - tdelay])
            x[p, i + 1] = x[p, i] + dt * ((1 / τd[p]) * (1 - x[p, i]) - r[p, i + 1] * x[p, i] * spd[i - tdelay])
            Is[p, i + 1] = Is[p, i] + dt * ((-Is[p, i] / τsI) + A[p] * r[p, i + 1] * x[p, i] * spd[i - tdelay])
        end
    end

    return vec(sum(Is, dims=1))
end

# Atividade de fundo poissoniana:
w_ps = 1.0
I_ps = zeros(6, 2, nSim)

if w_ps != 0
    W_ps = w_ps * rand(6, 2)

    for L in 1:6
        fr = 20 + 2 * randn()  # Frequência de disparo poissoniana
        spikess = poisson_spike_gen(fr, T / 1000, 1, dt / 1000)
        tps = findall(x -> x == 1, spikess[1, :]) 

        if !isempty(tps)
            I_ps[L, 1, :] .= W_ps[L, 1] * TMsynE(tps, nSim, τ_syn, dt)
            I_ps[L, 2, :] .= W_ps[L, 2] * TMsynI(tps, nSim, τ_syn, dt)
        end
    end
end

