function dbs_delta(fdbs, dbs_duration, dev, nSim, Fs, cut)
    Tdbs = 1 / fdbs  # Time interval between DBS pulses in seconds
    N_steps = Fs * Tdbs  # Number of steps between DBS pulses (Fs = 1000/dt Hz)
    dbs = 1:round(Int, N_steps):dbs_duration  # Indices of DBS pulses

    I_dbs_full = zeros(dbs_duration)
    I_dbs_full[dbs] .= 1  # Assigning Dirac delta pulses

    if dev == 1
        return I_dbs_full
    else
        return vcat(
            zeros(Int32(cut)),  # Initial zero padding
            zeros(Int32((nSim - cut) ÷ dev)),  # First inactive segment
            I_dbs_full,  # DBS active segment
            zeros(Int32((nSim - cut) ÷ dev))  # Final inactive segment
        )
    end
end

function TMsynE_dbs(dbs, nSim, tdelay, dt=0.1)
    # Synaptic parameters for 3 types: Facilitating, Depressing, Pseudo-linear
    τf = [670, 17, 326]    # Facilitation time constants (ms)
    τd = [138, 671, 329]   # Depression time constants (ms)
    U = [0.09, 0.5, 0.29]    # Utilization factor (release probability)
    A = [0.20, 0.63, 0.17]   # Synaptic efficacy
    τsE = 3               # Synaptic current decay time constant (ms)

    # Initialize state variables
    r = zeros(3, nSim)  # Recovered neurotransmitter resources
    x = ones(3, nSim)   # Available neurotransmitter resources
    Is = zeros(3, nSim) # Synaptic current contributions

    # Iterate over synapse types (F, D, P)
    for p in 1:3
        for i in (tdelay + 1):(nSim - 1)
            r[p, i + 1] = r[p, i] + dt * (-r[p, i] / τf[p] + U[p] * (1 - r[p, i]) * dbs[i - tdelay])
            x[p, i + 1] = x[p, i] + dt * ((1 / τd[p]) * (1 - x[p, i]) - r[p, i + 1] * x[p, i] * dbs[i - tdelay])
            Is[p, i + 1] = Is[p, i] + dt * ((-Is[p, i] / τsE) + A[p] * r[p, i + 1] * x[p, i] * dbs[i - tdelay])
        end
    end

    # Total postsynaptic current due to DBS
    I_dbs = sum(Is, dims=1)[1, :]

    return I_dbs
end

# Set the % of affected neurons in D layer by DBS
nh = 0.1     
# Impact of DBS on the other cortical structures via D PNs axons
fidCI = abs(1 * fidD)   # Synaptic fidelity for DBS carriers to invade CIs
fidM = abs(0 * fidD)    # Synaptic fidelity for DBS carriers to invade layer M
fidS = abs(1 * fidD)    # Synaptic fidelity for DBS carriers to invade layer S
fidR = abs(1 * fidD)    # Synaptic fidelity for DBS carriers to invade layer TCR
fidN = abs(1 * fidD)    # Synaptic fidelity for DBS carriers to invade layer TRN

nCI, nS, nR, nN = 1, 1, 1, 1
n_conn_CI = nCI * nh * nINs     # Percentage of CI neurons with synaptic contact from hyperdirect axon arbors
n_conn_S = nS * nh * nEs        # Percentage of S neurons with synaptic contact from hyperdirect axon arbors
n_conn_M = 0 * nh * nEm         # Percentage of M neurons with synaptic contact from hyperdirect axon arbors
n_conn_R = nR * nh * nErel      # Percentage of R neurons with synaptic contact from hyperdirect axon arbors
n_conn_N = nN * nh * nIret      # Percentage of N neurons with synaptic contact from hyperdirect axon arbors

n_hyp = nEd * nh  # Number of hyperdirect neurons

# Initialize DBS current
I_dbs = zeros(2, nSim)
if fidD != 0
    dev = 3
    fdbs = 130

    dbs_duration = (nSim - chop_till) ÷ dev  # in seconds
    I_dbs_pre = dbs_delta(fdbs, Int32(dbs_duration), dev, nSim, Fs, chop_till)

    # Postsynaptic DBS pulses (intra-axonal)
    I_dbs_post = TMsynE_dbs(I_dbs_pre, nSim, τ_syn, dt)

    I_dbs[1, :] .= I_dbs_pre
    I_dbs[2, :] .= I_dbs_post
end