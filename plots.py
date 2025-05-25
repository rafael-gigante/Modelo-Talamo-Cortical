import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import scipy.signal

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300
})

def plot_coupling_matrix(ZZ, Zx, condition="Normal"):
    dd = "plots"
    os.makedirs(dd, exist_ok=True)

    labels = ['S', 'M', 'P', 'INs', 'Rel', 'Ret']

    for matrix, suffix in zip([ZZ, Zx], ["", "_Normalized"]):
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.heatmap(matrix, annot=True, cmap='bwr', vmin=-1, vmax=1,
                    xticklabels=labels, yticklabels=labels)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)
        #plt.title(f"TCM Coupling Matrix ({condition}{suffix})", pad=20)
        plt.tight_layout()
        figdest = os.path.join(dd, f"TCM_Coupling_Matrix{suffix}_{condition}.png")
        plt.savefig(figdest, bbox_inches="tight")
        plt.close()

def dbs_plot(I_dbs, filename, xlim=None):
    plt.figure(figsize=(8, 4))
    t = np.arange(0, len(I_dbs)) * 0.1
    plt.plot(t, I_dbs, linewidth=1, color='k')
    if xlim:
        plt.xlim(xlim)
    plt.xlabel("Tempo (ms)")
    plt.ylabel("Amplitude da DBS")
    #plt.title("Sinal de Estimulação")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename + ".png", bbox_inches="tight")
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

    fig, ax = plt.subplots(figsize=(12, 6))
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


    ax.set_xlabel("Tempo (ms)")
    ax.set_ylabel("Neurônio")
    ax.set_ylim(0,n_tot)
    dd = "plots"
    if fidD == 0:
        #ax.set_title("Raster Plot - Condição PD")
        figdest = os.path.join(dd, f"raster_plot_PD")
    else:
        #ax.set_title("Raster Plot - Condição PD + DBS")
        figdest = os.path.join(dd, f"raster_plot_PD_DBS")
    plt.tight_layout()
    plt.savefig(figdest + ".png")
    plt.close()

def plot_PSD(f, S):
    plt.figure(figsize=(8, 4))
    plt.plot(f, S, c='k', linewidth=1)
    plt.axvspan(13, 30, color='grey', alpha=0.3, label='Banda Beta')
    plt.xlim([0, 100])
    plt.xticks(np.arange(0, 101, 10))
    plt.xlabel('Frequência (Hz)')
    plt.ylabel(r'PSD [$mV^2$/Hz]')
    #plt.title('Densidade Espectral de Potência (PSD)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig('plots/PSD.png', bbox_inches="tight")
    plt.close()

def plot_LFP(lfp, title):
    new_time = np.arange(len(lfp)) * 0.1
    plt.figure(figsize=(8, 4))
    plt.plot(new_time, lfp, c='k', linewidth=1)
    #plt.title(title)
    plt.xlabel('Tempo (ms)')
    plt.ylabel('LFP')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'plots/{title}.png', bbox_inches="tight")
    plt.close()

def PSD(signal, fs):
    nperseg = 10240  # 10 * 1024 para resolução alta
    f, S = scipy.signal.welch(signal, fs=fs, nperseg=nperseg, noverlap=nperseg//2, window='hann')
    return f, S

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = scipy.signal.butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = scipy.signal.filtfilt(b, a, data)  # filtfilt evita distorção de fase
    return y

# Carregamento dos dados
ZZ_diff = np.loadtxt("data/ZZ_diff.csv", delimiter=",")
Zx_diff = np.loadtxt("data/Zx_diff.csv", delimiter=",")
plot_coupling_matrix(ZZ_diff, Zx_diff, 'Difference')

I_dbs_pre = np.loadtxt("data/I_dbs_pre.csv", delimiter=",")
I_dbs_post = np.loadtxt("data/I_dbs_post.csv", delimiter=",")
dbs_plot(I_dbs_pre, "plots/DBS_pre")
dbs_plot(I_dbs_post, "plots/DBS_post")

dt = 0.1
nEs, nEm, nEd, nErel, nINs, nIret = 100, 100, 100, 100, 100, 40
n_tot = nEs + nEm + nEd + nINs + nErel + nIret

vEs = np.loadtxt("data/vEs.csv", delimiter=",")
vEm = np.loadtxt("data/vEm.csv", delimiter=",")
vEd = np.loadtxt("data/vEd.csv", delimiter=",")
vIs = np.loadtxt("data/vIs.csv", delimiter=",")
vErel = np.loadtxt("data/vErel.csv", delimiter=",")
vIret = np.loadtxt("data/vIret.csv", delimiter=",")

cEs = np.loadtxt("data/cEs.csv", delimiter=",")
cEm = np.loadtxt("data/cEm.csv", delimiter=",")
cEd = np.loadtxt("data/cEd.csv", delimiter=",")
cIs = np.loadtxt("data/cIs.csv", delimiter=",")
cErel = np.loadtxt("data/cErel.csv", delimiter=",")
cIret = np.loadtxt("data/cIret.csv", delimiter=",")

nSim = len(vEs[0])
fidD = 0  # 0 para PD, 1 para PD + DBS

raster_plot(nEs, nEm, nEd, nINs, nErel, nIret, n_tot,
            vEs, vEm, vEd, vIs, vIret, vErel,
            cIret, cErel, cIs, cEd, cEm, cEs,
            nSim, dt, fidD)

LFP = np.loadtxt("data/LFP.csv", delimiter=",")
lowcut = 13
highcut = 30
fs = 1000 / dt

LFP = butter_bandpass_filter(LFP, lowcut, highcut, fs, order=3)
plot_LFP(LFP, "LFP")

frequencies, psd = PSD(LFP, fs=fs)
plot_PSD(frequencies, psd)