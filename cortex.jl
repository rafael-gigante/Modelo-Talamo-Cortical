module Cortex
export S_layer, M_layer, D_layer, ctx_INs

using LinearAlgebra
include("TM_synapses.jl")
using .TMsynapses

function S_layer(a, b, c, d, n, v, u, r, x, Is,
                EPSC, EPSCd, EPSCm, EPSC_rel, IPSC_in, IPSC_ret,
                W_EEss, W_EEsd, W_EEsm, W_Erel, W_EI, W_Eret,
                I_psE, I_psI, kisi, zeta, Idc, Idbs, n_affected, dt)

    vp = 30
    sp = zeros(3)
    Ise = zeros(n)

    for k in 1:n
        if (k >= 1 && k < n_affected)
            Idbss = Idbs
        else
            Idbss = 0
        end

        v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
                      (W_EEss[k] * EPSC / n) +
                      (W_EEsd[k] * EPSCd / n) + (W_EEsm[k] * EPSCm / n) + (W_Erel[k] * EPSC_rel / n) +
                      (W_EI[k] * IPSC_in / n) + (W_Eret[k] * IPSC_ret / n) +
                      I_psE - I_psI + Idbss + kisi[k])
        
        u[k] += dt * a[k] * (b[k] * v[k] - u[k])
        
        if v[k] >= vp + zeta[k]
            v[k] = c[k]
            u[k] += d[k]
            sp .= 1
        end

        spikeE = sp
        r, x, Is, Ipost = TMsynE_inst(r, x, Is, spikeE)
        Ise[k] = Ipost

        sp .= 0
    end

    EPSC = sum(Ise)
    return v, u, r, x, Is, EPSC
end

function M_layer(a,b,c,d,n,v,u,r,x,Is,
                EPSC,EPSCs,EPSCd,EPSC_rel,IPSC_in,IPSC_ret,
                W_EEmm,W_EEms,W_EEmd,W_Erel,W_EI,W_Eret,
                I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
    vp = 30
    sp = zeros(3)
    Ise = zeros(n)

    for k in 1:n
        if (k >= 1 && k < n_affected)
            Idbss = Idbs
        else
            Idbss = 0
        end

        v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
                (W_EEmm[k] * EPSC / n) +  # Self feedback
                (W_EEms[k] * EPSCs / n) + (W_EEmd[k] * EPSCd / n) + (W_Erel[k] * EPSC_rel / n) + # Excitatory inputs
                (W_EI[k] * IPSC_in / n) + (W_Eret[k] * IPSC_ret / n) + # Inhibitory inputs
                I_psE - I_psI + Idbss + kisi[k])
        
        u[k] += dt * a[k] * (b[k] * v[k] - u[k])
        
        if v[k] >= vp + zeta[k]
            v[k] = c[k]
            u[k] += d[k]
            sp .= 1
        end

        spikeE = sp
        r, x, Is, Ipost = TMsynE_inst(r, x, Is, spikeE)
        Ise[k] = Ipost

        sp .= 0
    end

    EPSC = sum(Ise)
    return v, u, r, x, Is, EPSC
end

function D_layer(a, b, c, d, n, n_hyp, v, u, r, x, Is, rf, xf, Isf,
    EPSC, EPSCs, EPSCm, EPSC_rel, IPSC_in, IPSC_ret,
    W_EEdd, W_EEds, W_EEdm, W_Erel, W_EI, W_Eret,
    I_psE, I_psI, kisi, zeta, Idc, Idbs, dt)

    vp = 30.0
    sp = zeros(3)
    Ise = zeros(n)
    IseF = zeros(n)

    for k in 1:n 
        if n_hyp == 0
            Idbss = 0
        else
            if k >= n_hyp
                Idbss = Idbs[2] 
            else
                Idbss = Idbs[1]
            end
        end

    v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
            (W_EEdd[k] * EPSC / n) +  # Self feedback
            (W_EEds[k] * EPSCs / n) + (W_EEdm[k] * EPSCm / n) + (W_Erel[k] * EPSC_rel / n) +  # Excitatory inputs
            (W_EI[k] * IPSC_in / n) + (W_Eret[k] * IPSC_ret / n) + # Inhibitory inputs
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

    rsF, xsF, IsynF, IpostF = TMsynE_inst_F(rr, xx, Iss, spikeE)
    rf, xf, Isf = rsF, xsF, IsynF
    IseF[k] = IpostF

    sp .= 0  # Reset spike array
    end

    EPSC = sum(Ise)
    EPSCf = sum(IseF)

return v, u, r, x, Is, EPSC, rf, xf, Isf, EPSCf
end

function ctx_INs(a,b,c,d,n,v,u,r,x,Is,
                IPSC,EPSCs,EPSCm,EPSCd,EPSC_rel,IPSC_ret,
                W_IEd,W_IEs,W_IEm,W_Iret,W_II,W_Irel,
                I_psE,I_psI,kisi,zeta,Idc,Idbs,n_affected,dt)
    vp = 30
    sp = zeros(3)
    Ise = zeros(n)
            
    for k in 1:n
        if (k >= 1 && k < n_affected)
            Idbss = Idbs
        else
            Idbss = 0
        end

        v[k] += dt * (0.04 * v[k]^2 + 5 * v[k] - u[k] + 140 + Idc[k] +
                    (W_II[k] * IPSC / n) +  # Self feedback
                    (W_IEd[k] * EPSCd / n) + (W_IEs[k] * EPSCs / n) + (W_IEm[k] * EPSCm / n) + (W_Irel[k] * EPSC_rel / n) +  # Excitatory inputs
                    (W_Iret[k] * IPSC_ret / n) +  # Inhibitory inputs
                    I_psE - I_psI + Idbss + kisi[k])
        u[k] += dt * a[k] * (b[k] * v[k] - u[k])
        
        if v[k] >= vp + zeta[k]
            v[k] = c[k]
            u[k] += d[k]
            sp .= 1
        end
            
        spikeE = sp
        r, x, Is, Ipost = TMsynE_inst(r, x, Is, spikeE)
        Ise[k] = Ipost
            
        sp .= 0
    end
            
    EPSC = sum(Ise)
    return v, u, r, x, Is, EPSC
end

end