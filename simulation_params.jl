# SIMULATION TIME AND STEP, Time delays       
T = (sim_time + 1) * 1000 #total simulation time in ms with and additional second
dt = 0.1                  #time step in ms
Fs = Int32(1000 / dt)     #sampling frequency in Hz
nSim = round(Int, T / dt) #number of simulation steps
chop_till = 1 * Fs        #chop the first second of the simulation
τ_L = 8                   #time delay between different strucutures within the cortex and thalamus in ms
τ_wL = 1                  #time delay within the same structure in ms     
τ_TC = 15                 #time delay from thalamus to cortex in ms
τ_CT = 20                 #time delay from cortex to thalamus in ms
τ_syn = 1                 #synaptic transmission delay for all synapses

# Define the time vector
if τ_TC >= τ_CT
    tVec = (τ_TC + τ_syn + 1):nSim-1
else
    tVec = (τ_CT + τ_syn + 1):nSim-1
end