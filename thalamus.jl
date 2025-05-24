module Thalamus

export thm_rel, thm_ret

using LinearAlgebra
include("TM_synapses.jl")
using .TMsynapses

function thm_rel(a, b, c, d, n, v, u, r, x, Is, rd, xd, Isd,
                 EPSC, EPSCs, EPSCm, EPSCd, IPSC_in, IPSC_ret,
                 W_EE, W_EErs, W_EErm, W_EErd, W_EI_IN, W_EI_ret,
                 I_psE, I_psI, kisi, zeta, Idc, Idbs, n_affected, dt)
    
    vp = 30
    sp = zeros(3)
    Ise = zeros(n)
    IseD = zeros(n)
    
    for k in 1:n
        if (k >= 1 && k < n_affected)
            Idbss = Idbs
        else
            Idbss = 0
        end
        
        v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
                      (W_EE[k] * EPSC / n) +  # Self feedback
                      (W_EErs[k] * EPSCs / n) + (W_EErm[k] * EPSCm / n) + (W_EErd[k] * EPSCd / n) +  # Excitatory inputs
                      (W_EI_ret[k] * IPSC_ret / n) + (W_EI_IN[k] * IPSC_in / n) +  # Inhibitory inputs
                      I_psE - I_psI + Idbss + kisi[k])
        u[k] += dt * a[k] * (b[k] * v[k] - u[k])
        
        if v[k] >= vp + zeta[k]
            v[k] = c[k]
            u[k] += d[k]
            sp .= 1
        end
        
        spikeE = sp
        rr, xx, Iss = r, x, Is

        rs, xs, Isyn, Ipost = TMsynE_inst(r, x, Is, spikeE)
        r, x, Is = rs, xs, Isyn
        Ise[k] = Ipost
        
        rsD, xsD, IsynD, IpostD = TMsynE_inst_D(rr, xx, Iss, spikeE)
        rd, xd, Isd = rsD, xsD, IsynD
        IseD[k] = IpostD
        
        sp .= 0
    end
    
    EPSC = sum(Ise)
    EPSC_d = sum(IseD)
    return v, u, r, x, Is, EPSC, rd, xd, Isd, EPSC_d
end

function thm_ret(a, b, c, d, n, v, u, r, x, Is,
                 IPSC, EPSCs, EPSCm, EPSCd, IPSC_in, EPSC_rel,
                 W_II, W_IErs, W_IErm, W_IErd, W_II_IN, W_IE_rel,
                 I_psE, I_psI, kisi, zeta, Idc, Idbs, n_affected, dt)
    
    vp = 30
    sp = zeros(3)
    Isi = zeros(n)
    
    for k in 1:n
        if (k >= 1 && k < n_affected)
            Idbss = Idbs
        else
            Idbss = 0
        end
        
        v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
                      (W_II[k] * IPSC / n) +  # Self feedback
                      (W_IErs[k] * EPSCs / n) + (W_IErm[k] * EPSCm / n) + (W_IErd[k] * EPSCd / n) + (W_IE_rel[k] * EPSC_rel / n) +  # Excitatory inputs
                      (W_II_IN[k] * IPSC_in / n) +  # Inhibitory inputs
                      I_psE - I_psI + Idbss + kisi[k])
        u[k] += dt * a[k] * (b[k] * v[k] - u[k])
        
        if v[k] >= vp + zeta[k]
            v[k] = c[k]
            u[k] += d[k]
            sp .= 1
        end
        
        spikeI = sp
        rs, xs, Isyn, Ipost = TMsynI_inst(r, x, Is, spikeI)
        r, x, Is = rs, xs, Isyn
        Isi[k] = Ipost
        
        sp .= 0
    end
    
    IPSC = sum(Isi)
    return v, u, r, x, Is, IPSC
end
end