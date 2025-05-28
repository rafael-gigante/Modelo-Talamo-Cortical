import matplotlib.pyplot as plt
import seaborn as sns
import os
import csv
import numpy as np
import scipy.signal
from scipy.stats import ks_2samp
from scipy.integrate import simpson

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
    windowlen = 4 * fs 
    f, S = scipy.signal.welch(signal,fs=fs,nperseg=windowlen)
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

with open("psd_analysis.csv", mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['normal', "PD", "PD + DBS"])

for i in range(1,11):
    LFP_dbs = np.loadtxt(f"dbs_n/trial_{i}/LFP.csv", delimiter=",")
    LFP_norm = np.loadtxt(f"normal_n/trial_{i}/LFP.csv", delimiter=",")
    LFP_PD = np.loadtxt(f"dp_n/trial_{i}/LFP.csv", delimiter=",")

    lowcut = 13
    highcut = 30
    fs = 1000/0.1

    LFP_norm = butter_bandpass_filter(LFP_norm, lowcut, highcut, fs, order=3)
    LFP_PD = butter_bandpass_filter(LFP_PD, lowcut, highcut, fs, order=3)
    LFP_dbs = butter_bandpass_filter(LFP_dbs, lowcut, highcut, fs, order=3)

    intervalos = [
    [[0, 10000], [50000, 60000], [100000, 110000]],  # Primeiro bloco
    [[10000, 20000], [60000, 70000], [110000, 120000]],   # Segundo bloco
    [[20000, 30000], [70000, 80000], [120000, 130000]],    # Terceiro bloco
    [[30000, 40000], [80000, 90000], [130000, 140000]],    # Quarto bloco
    [[40000, 50000], [90000, 100000], [140000, 150000]]    # Quinto bloco
    ]


    for j in range(5):

        LFP_norm_aux = LFP_norm[intervalos[j][0][0]:intervalos[j][0][1]]
        LFP_PD_aux = LFP_PD[intervalos[j][1][0]:intervalos[j][1][1]]
        LFP_dbs_aux = LFP_dbs[intervalos[j][2][0]:intervalos[j][2][1]]

        frequencies_norm, psd_norm = PSD(LFP_norm_aux, fs=fs)
        frequencies_PD, psd_PD = PSD(LFP_PD_aux, fs=fs)
        frequencies_dbs, psd_dbs = PSD(LFP_dbs_aux, fs=fs)

        p_norm = simpson(y=psd_norm, x=frequencies_norm)
        p_PD = simpson(y=psd_PD, x=frequencies_PD)
        p_dbs = simpson(y=psd_dbs, x=frequencies_dbs)

        with open("psd_analysis.csv", mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([p_norm, p_PD, p_dbs])   


plt.figure(figsize=(8, 4))
plt.plot(frequencies_norm, psd_norm, c='b', linewidth=1, label='Normal')
plt.plot(frequencies_PD, psd_PD, c='r', linewidth=1, label='DP')
plt.plot(frequencies_dbs, psd_dbs, c='g', linewidth=1, label='DP + ECP (130 Hz)')
plt.axvspan(13, 30, color='grey', alpha=0.3, label=r'Banda $\beta$ (13-30 Hz)')
plt.xlim([0, 100])
plt.xticks(np.arange(0, 101, 10))
plt.xlabel('Frequência (Hz)')
plt.ylabel(r'PSD [$mV^2$/Hz]')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('plots/PSD.png', bbox_inches="tight")
plt.close()