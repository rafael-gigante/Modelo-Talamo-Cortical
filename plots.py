import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import matplotlib.mlab as mlab
import scipy.signal 

def plot_coupling_matrix(ZZ, Zx, condition="Normal"):
    dd = "Modelo TC - julia\plots"
    if not os.path.exists(dd):
        os.makedirs(dd)

    # Plot and save the unnormalized matrix
    plt.figure()
    ax = sns.heatmap(ZZ, annot=True, cmap='bwr', vmin=-1, vmax=1)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.xticks(np.arange(0.5, len(ZZ), 1), labels=['S', 'M', 'D', 'CI', 'TCR', 'TRN'])
    plt.yticks(np.arange(0.5, len(ZZ), 1), labels=['S', 'M', 'D', 'CI', 'TCR', 'TRN'], rotation='horizontal')
    plt.title(f"TCM Coupling Matrix ({condition})")
    figdest = os.path.join(dd, f"TCM_Coupling_Matrix_{condition}")
    plt.savefig(figdest + ".png")
    plt.close()

    # Plot and save the normalized matrix
    plt.figure()
    ax = sns.heatmap(Zx, annot=True, cmap='bwr', vmin=-1, vmax=1)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.xticks(np.arange(0.5, len(Zx), 1), labels=['S', 'M', 'D', 'CI', 'TCR', 'TRN'])
    plt.yticks(np.arange(0.5, len(Zx), 1), labels=['S', 'M', 'D', 'CI', 'TCR', 'TRN'], rotation='horizontal')
    plt.title(f"TCM Coupling Matrix ({condition})")
    figdest = os.path.join(dd, f"TCM_Coupling_Matrix_Normalized_{condition}")
    plt.savefig(figdest + ".png")
    plt.close()

def dbs_plot(I_dbs, filename, xlim=None):
    plt.figure()
    t = np.arange(0, len(I_dbs))
    plt.plot(t*0.1, I_dbs, linewidth=1, color='k')
    if xlim:
        plt.xlim(xlim)
    plt.xlabel("Time (ms)")
    plt.ylabel("DBS Amplitude")
    plt.title(filename.split("\\")[-1])  # Extracting file name for title
    plt.grid(True)
    # Save as .png (can change to .svg or .jpeg)
    plt.savefig(filename + ".png", dpi=300, bbox_inches="tight")
    plt.close()

def raster_plot(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,vEs,vEm,vEd,vIs,vRet,vRel,cIret,cErel,cIs,cEd,cEm,cEs,nSim,dt, fidD):
    def get_spikes(v, c, offset, dt=0.1):
        for neuron in range(len(v)):
            fired = np.where(v[neuron] == c[neuron])[0]
            firings.extend([[t * dt, neuron + offset] for t in fired])
    
    firings = []
    offsets = np.cumsum([0, nIret, nErel, nINs, nEd, nEm])
    voltages = [vRet, vRel, vIs, vEd, vEm, vEs]
    conductances = [cIret, cErel, cIs, cEd, cEm, cEs]
    for i in range(len(voltages)):
        get_spikes(voltages[i], conductances[i], offsets[i])

    firings = np.array(firings)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(firings[:, 0], firings[:, 1], s=1, color='black')

    # Definir separações visuais para as estruturas
    magenta_lines = [nIret+nErel+nINs/2,  nIret + nErel + nINs + (0.7*nEd), nIret + nErel + nINs + nEd + nEm + nEs/2]

    structure_divison = [0, nIret, nIret + nErel, nIret + nErel + nINs, nIret + nErel + nINs + nEd, 
                nIret + nErel + nINs + nEd + nEm, n_tot]

    structures_text_pos = [nIret/2, nIret + nErel/2, nIret + nErel + nINs/2, nIret + nErel + nINs + nEd/2, 
                nIret + nErel + nINs + nEd + nEm/2, n_tot - nEs/2]

    neuron_types_division = [0, nIret, nIret + nErel, nIret + nErel + (nINs/2), nIret + nErel + nINs, 
                nIret + nErel + nINs + (0.7*nEd), nIret + nErel + nINs + nEd, 
                nIret + nErel + nINs + nEd + nEm, n_tot-(nEs/2), n_tot]

    neuron_types_text_pos = [nIret/2, nIret + nErel/2, nIret + nErel + (nINs/2)/2, nIret + nErel + nINs - (nINs/2)/2, 
                nIret + nErel + nINs + (0.7*nEd)/2, nIret + nErel + nINs + nEd - (0.3*nEd)/2, 
                nIret + nErel + nINs + nEd + nEm/2, nIret + nErel + nINs + nEd + nEm + (nEs/2)/2, n_tot - (nEs/2)/2]    


    color1 = ['orange', 'cyan', 'green', 'magenta', 'blue', 'red']
    color2 = ['gold', 'blue', 'cyan', 'magenta', 'green', 'red', 'green', 'green', 'red']
    structures = ['TRN', 'TCR', 'CI', 'D', 'M', 'S']
    neuron_types = ['TR', 'TC', 'FS', 'LTS', 'RS', 'IB', 'RS', 'RS', 'IB']

    t_tot = nSim*dt
    #for i in magenta_lines:
    #    ax.axhline(i, xmin=(t_tot*0.05)/(t_tot+((t_tot*0.1))), xmax=((t_tot)+((t_tot*0.05)))/(t_tot+(t_tot*0.1)), color='magenta', linestyle='--', linewidth=2)

    for i in range(len(structures)):
        ax.text(-(0.1*t_tot*0.25), structures_text_pos[i], structures[i], ha='center', va='center', fontsize=12, fontstyle='italic', fontweight='bold', color='k')

    for i in range(len(neuron_types)):
        ax.text((t_tot)+(0.25*0.1*t_tot), neuron_types_text_pos[i], neuron_types[i], ha='center', va='center', fontsize=12, fontstyle='italic', fontweight='bold', color='k')


    for i in range(len(structure_divison)-1):
        ax.axvline(x=0, ymin=structure_divison[i]/n_tot, ymax=structure_divison[i+1]/n_tot, color=color1[i], linestyle='-', linewidth=2)
        #ax.axhline(structure_divison[i+1], xmin=(t_tot*0.05)/(t_tot+((t_tot*0.1))), xmax=((t_tot)+((t_tot*0.05)))/(t_tot+(t_tot*0.1)),color='k', linestyle='-', linewidth=2)

    for i in range(len(neuron_types_division)-1):
        ax.axvline(x=(t_tot), ymin=neuron_types_division[i]/n_tot, ymax=neuron_types_division[i+1]/n_tot, color=color2[i], linestyle='-', linewidth=2)


    ax.set_xlabel("time (ms)")
    ax.set_ylabel("neuron #")
    ax.set_ylim(0,n_tot)
    dd = "Modelo TC - julia\plots"
    if fidD == 0:
        ax.set_title("Neuronal Spiking Raster Plot - PD condition")
        figdest = os.path.join(dd, f"raster_plot_PD")
    else:
        ax.set_title("Neuronal Spiking Raster Plot - PD + DBS condition")
        figdest = os.path.join(dd, f"raster_plot_PD_DBS")
    plt.tight_layout()
    plt.savefig(figdest + ".png")
    plt.close()

    

ZZ_diff = np.loadtxt("Modelo TC - julia/data/ZZ_diff.csv", delimiter=",")
Zx_diff = np.loadtxt("Modelo TC - julia/data/Zx_diff.csv", delimiter=",")
plot_coupling_matrix(ZZ_diff, Zx_diff, 'Difference')

I_dbs_pre = np.loadtxt("Modelo TC - julia/data/I_dbs_pre.csv", delimiter=",")
I_dbs_post = np.loadtxt("Modelo TC - julia/data/I_dbs_post.csv", delimiter=",")

dt = 0.1

nEs, nEm, nEd, nErel, nINs, nIret = 100, 100, 100, 100, 100, 40
n_tot = nEs + nEm + nEm + nINs + nErel + nIret
vEs = np.loadtxt("Modelo TC - julia/data/vEs.csv", delimiter=",")
vEm = np.loadtxt("Modelo TC - julia/data/vEm.csv", delimiter=",")
vEd = np.loadtxt("Modelo TC - julia/data/vEd.csv", delimiter=",")
vIs = np.loadtxt("Modelo TC - julia/data/vIs.csv", delimiter=",")
vErel = np.loadtxt("Modelo TC - julia/data/vErel.csv", delimiter=",")
vIret = np.loadtxt("Modelo TC - julia/data/vIret.csv", delimiter=",")
cEs = np.loadtxt("Modelo TC - julia/data/cEs.csv", delimiter=",")
cEm = np.loadtxt("Modelo TC - julia/data/cEm.csv", delimiter=",")
cEd = np.loadtxt("Modelo TC - julia/data/cEd.csv", delimiter=",")
cIs = np.loadtxt("Modelo TC - julia/data/cIs.csv", delimiter=",")
cErel = np.loadtxt("Modelo TC - julia/data/cErel.csv", delimiter=",")
cIret = np.loadtxt("Modelo TC - julia/data/cIret.csv", delimiter=",")
nSim = len(vEs[0])
fidD = 0 # 0 for PD condition, 1 for PD + DBS condition
raster_plot(nEs,nEm,nEd,nINs,nErel,nIret,n_tot,vEs,vEm,vEd,vIs,vIret,vErel,cIret,cErel,cIs,cEd,cEm,cEs,nSim,dt,fidD)

def PSD(signal, fs):
    (f, S) = scipy.signal.welch(signal, fs, nperseg=10*1024)
    
    return f, S

def plot_PSD(f, S):
    x_arr = np.arange(0, 101, 10)
    
    plt.figure()
    plt.plot(f, S, c='k', linewidth=1)
    plt.axvspan(13, 30, color='grey', alpha=0.3)
    plt.xlim([0, 100])
    plt.xticks(x_arr)
    plt.xlabel('frequency (Hz)')
    plt.ylabel(r'PSD [$mV^2/Hz$]')
    plt.title('PSD')
    plt.savefig('Modelo TC - julia\plots\PSD.png')
    plt.close()
    
def plot_LFP(lfp, title):
    new_time= np.transpose(np.arange(len(lfp)))
    
    plt.figure()
    
    plt.title(title)

    plt.plot(new_time*0.1, lfp, c='k', linewidth=1)
    
    # Set the x-axis label
    plt.xlabel('Time (ms)')
    plt.ylabel('LFP')
    
    plt.savefig(f'Modelo TC - julia\plots\{title}.png')
    plt.close()


LFP = np.loadtxt("Modelo TC - julia/data/LFP.csv", delimiter=",")
plot_LFP(LFP, "LFP")

# Compute Power Spectral Density using Welch's method
frequencies, psd = PSD(LFP, fs=1000/dt)
plot_PSD(frequencies, psd)