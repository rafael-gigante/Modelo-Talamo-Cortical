# Parâmetros da simulação:
T = (sim_time + 1) * 1000 # tempo total de simulação em ms com 1 segundo adicional
dt = 0.1                  # passo de tempo em ms
Fs = Int32(1000 / dt)     # frequência de amostragem em Hz
nSim = round(Int, T / dt) # número de passos de simulação
chop_till = 1 * Fs        # descartar o primeiro segundo da simulação
τ_L = Int32(8 / dt)                   # atraso de tempo entre diferentes estruturas dentro do córtex e tálamo em ms
τ_wL = Int32(1 / dt)                  # atraso de tempo dentro da mesma estrutura em ms
τ_TC = Int32(15 / dt)                 # atraso de tempo do tálamo para o córtex em ms
τ_CT = Int32(20 / dt)                 # atraso de tempo do córtex para o tálamo em ms
τ_syn = Int32(1 / dt)                 # atraso de transmissão sináptica para todas as sinapses

# Define o vetor com os passos de tempo da simulação
if τ_TC >= τ_CT
    tVec = Int32(τ_TC + τ_syn + 1/dt):nSim-1
else
    tVec = Int32(τ_CT + τ_syn + 1/dt):nSim-1
end