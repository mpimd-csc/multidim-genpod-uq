import numpy as np
import matplotlib.pyplot as plt

import gen_pod_uq.mc_pce_utils as mpu

from kmg_mtrx_helpers import get_evay

smpls = 1  # whether to average the rb/rand bases

def comp_cdf(ysoltens, nua=0., nub=1., dst='beta-2-5', nsample=200000,
             pcedim=5):
    
    abscissae, _, compexpv, _ = mpu.\
        setup_pce(distribution=dst,
                  distrpars=dict(a=nua, b=nub),
                  pcedim=pcedim, uncdims=4)
    
    expv = compexpv(ysoltens)
    evay = get_evay(ysoltens, abscissae)
    # vecy = ysoltens.reshape((1, -1))
    
    print(f'expv: {expv.item():.4e}')

    # yamn = evay([numean]*4)
    # print(f'y(amean): {yamn:.4e}')
    # yamn = evay([nua]*4)
    # print(f'y(amin): {yamn:.4e}')
    # yamn = evay([nub]*4)
    # print(f'y(amax): {yamn:.4e}')
    
    getsample = mpu.get_nu_sample(distribution=dst,
                                  uncdims=4, nulb=nua, nuub=nub)
    # randa = getsample(1)
    # yamn = evay(randa.flatten())
    # print(f'y(arnd): {yamn:.4e}')
    
    rndsa = getsample(nsample)
    smpllist = []
    for csmpl in rndsa:
        smpllist.append(evay(csmpl.flatten()))
    
    cpfvals = mpu.empirical_cdf(smpllist)
    srtdsmpllist = sorted(smpllist)

    return srtdsmpllist, cpfvals

if __name__ == '__main__':
    np.random.seed(1)

    nua, nub = 5e-4, 10e-4
    dst = 'uniform'
    dst = 'beta-2-5'
    Nndstr = f'N12nu{nua:.2e}--{nub:.2e}' + dst
    dataprfx = 'mh-data/cached-data/'  + Nndstr

    yts = dataprfx + '_pce5_ysoltns.npy'
    ysoltens = np.load(yts)
    xtrth, cdfxtrth = comp_cdf(ysoltens, pcedim=5, dst=dst, nua=nua, nub=nub)
    jmin, jmax = xtrth[0], xtrth[-1]

    yts = dataprfx + '_pce2_ysoltns.npy'
    ysoltens = np.load(yts)
    xpcetwo, cdfpcetwo = comp_cdf(ysoltens, pcedim=2, dst=dst, nua=nua, nub=nub)
    jmin, jmax = max(jmin, xpcetwo[0]), min(jmax, xpcetwo[-1])

    yts = dataprfx + '_pce5_pod8_bfpce2_run1of1_ysoltns.npy'
    ysoltens = np.load(yts)
    xpodpcef, cdfpodpcef = comp_cdf(ysoltens, pcedim=5, dst=dst, nua=nua, nub=nub)
    jmin, jmax = max(jmin, xpodpcef[0]), min(jmax, xpodpcef[-1])

    accytns = 0
    for kkk in range(smpls):
        cyts = dataprfx + '_pce5_pod8_bfrb_random16_runs10' + \
                f'_run{kkk+1}of10_ysoltns.npy'
        accytns += np.load(cyts)
    accytns = 1/smpls*accytns
    xrb, cdfrbx = comp_cdf(accytns, pcedim=5, dst=dst, nua=nua, nub=nub)
    jmin, jmax = max(jmin, xrb[0]), min(jmax, xrb[-1])

    accytns = 0
    for kkk in range(smpls):
        cyts = dataprfx + '_pce5_pod8_bfmc16_runs10' + \
                f'_run{kkk+1}of10_ysoltns.npy'
        accytns += np.load(cyts)
    accytns = 1/smpls*accytns
    xmc, cdfmcx = comp_cdf(accytns, pcedim=5, dst=dst, nua=nua, nub=nub)
    jmin, jmax = max(jmin, xmc[0]), min(jmax, xmc[-1])

    
    # plt.figure(1)
    # plt.plot(xtrth, cdfxtrth, label='FOM PCE[5]')
    # plt.plot(xpcetwo, cdfpcetwo, label='FOM PCE[2]')
    # plt.legend()

    iabsc = np.linspace(jmin, jmax, 1000)
    icdftrth = np.interp(iabsc, xtrth, cdfxtrth)
    icdfpcetwo = np.interp(iabsc, xpcetwo, cdfpcetwo)
    icdfpodpcef = np.interp(iabsc, xpodpcef, cdfpodpcef)
    icdfrb = np.interp(iabsc, xrb, cdfrbx)
    icdfmc = np.interp(iabsc, xmc, cdfmcx)

    # plt.figure(2)
    # plt.plot(iabsc, icdftrth, label='FOM PCE[5]')
    # plt.plot(iabsc, icdfpcetwo, label='FOM PCE[2]')
    # plt.legend()

    plt.figure(3)
    # plt.plot(iabsc, np.abs(icdfpcetwo-icdftrth), label='PCE[2]')
    plt.plot(iabsc, np.abs(icdfpodpcef-icdftrth), label='pcePOD16')
    plt.plot(iabsc, np.abs(icdfrb-icdftrth), label='wRB16')
    plt.plot(iabsc, np.abs(icdfmc-icdftrth), label='rndPOD16')
    plt.legend()

    plt.show()
