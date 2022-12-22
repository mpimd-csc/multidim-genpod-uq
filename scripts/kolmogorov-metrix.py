import numpy as np
import matplotlib.pyplot as plt

import gen_pod_uq.mc_pce_utils as mpu

from kmg_mtrx_helpers import get_evay


def comp_cdf(ysoltens, nua=0., nub=1., dst='beta-2-5', nsample=1000000,
             pcedim=5):
    
    abscissae, _, _, _ = mpu.\
        setup_pce(distribution=dst,
                  distrpars=dict(a=nua, b=nub),
                  pcedim=pcedim, uncdims=4)
    
    # expv = compexpv(ysoltens)
    evay = get_evay(ysoltens, abscissae)
    # vecy = ysoltens.reshape((1, -1))
    # print(f'expv: {expv.item():.4e}')

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

def compmaxdiff(xl, cdfxl, tx, tcdfx, intpoints=2000):
    smin, smax = tx[0], tx[-1]
    for x in xl:
        smin = max(smin, x[0])
        smax = min(smax, x[-1])
    intx = np.linspace(smin, smax, intpoints)
    itcdf = np.interp(x=intx, xp=tx, fp=tcdfx)

    diffl, maxl = [], []
    for kkk, cdfx in enumerate(cdfxl):
        icdf = np.interp(x=intx, xp=xl[kkk], fp=cdfx)
        dficdf = icdf-itcdf
        diffl.append([intx, dficdf])
        maxl.append(np.max(np.abs(dficdf)))

    return np.median(np.array(maxl)), maxl, diffl


if __name__ == '__main__':

    smpls = 10  # number of samples for the MC/wRB bases
    runs = 5  # how many runs --- since the sampling is also stochastic
    # np.random.seed(1)

    nua, nub = 5e-4, 10e-4
    dst = 'beta-2-5'
    dst = 'uniform'
    Nndstr = f'N12nu{nua:.2e}--{nub:.2e}' + dst
    dataprfx = 'mh-data/cached-data/'  + Nndstr

    for rrr in range(runs):
        yts = dataprfx + '_pce5_ysoltns.npy'
        ysoltens = np.load(yts)
        xtrth, cdfxtrth = comp_cdf(ysoltens, pcedim=5, dst=dst, nua=nua, nub=nub)
        # jmin, jmax = xtrth[0], xtrth[-1]

        yts = dataprfx + '_pce2_ysoltns.npy'
        ysoltens = np.load(yts)
        xpcetwo, cdfpcetwo = comp_cdf(ysoltens, pcedim=2, dst=dst, nua=nua, nub=nub)
        ppkmmed, _, ppxc \
            = compmaxdiff([xpcetwo], [cdfpcetwo], xtrth, cdfxtrth)
        print(f'Kolmometer: pce[2]: {ppkmmed:.5f}')
        # jmin, jmax = max(jmin, xpcetwo[0]), min(jmax, xpcetwo[-1])

        yts = dataprfx + '_pce5_pod8_bfpce2_run1of1_ysoltns.npy'
        ysoltens = np.load(yts)
        xpodpcef, cdfpodpcef = comp_cdf(ysoltens, pcedim=5, dst=dst, nua=nua, nub=nub)
        # jmin, jmax = max(jmin, xpodpcef[0]), min(jmax, xpodpcef[-1])
        ppkmmed, _, ppxc \
            = compmaxdiff([xpodpcef], [cdfpodpcef], xtrth, cdfxtrth)
        print(f'Kolmometer: pce-16: {ppkmmed:.5f}')

        # accytns = 0
        xrbl, rbcdfl = [], []
        for kkk in range(smpls):
            cyts = np.load(dataprfx + '_pce5_pod8_bfrb_random16_runs10' + \
                           f'_run{kkk+1}of10_ysoltns.npy')
            xrb, cdfrbx = comp_cdf(cyts, pcedim=5, dst=dst, nua=nua, nub=nub)
            xrbl.append(xrb)
            rbcdfl.append(cdfrbx)
            # accytns += cyts

        rbkmmed, rbkmerrs, rbxc = compmaxdiff(xrbl, rbcdfl, xtrth, cdfxtrth)
        print(f'Kolmometer: rb16: {rbkmmed:.5f} -- median out of {smpls}')

        xrblt, rbcdflt = [], []
        for kkk in range(smpls):
            cyts = np.load(dataprfx + '_pce5_pod8_bfrb_random32_runs10' + \
                           f'_run{kkk+1}of10_ysoltns.npy')
            xrbt, cdfrbxt = comp_cdf(cyts, pcedim=5, dst=dst, nua=nua, nub=nub)
            xrblt.append(xrbt)
            rbcdflt.append(cdfrbxt)
            # accytns += cyts

        rbkmmedt, _, _ = compmaxdiff(xrblt, rbcdflt, xtrth, cdfxtrth)
        print(f'Kolmometer: rb32: {rbkmmedt:.5f} -- median out of {smpls}')

        # accytns = 1/smpls*accytns
        # jmin, jmax = max(jmin, xrb[0]), min(jmax, xrb[-1])

        # accytns = 0
        xmcl, mccdfl = [], []
        for kkk in range(smpls):
            cyts = np.load(dataprfx + '_pce5_pod8_bfmc16_runs10' + \
                           f'_run{kkk+1}of10_ysoltns.npy')
            xmc, cdfmcx = comp_cdf(cyts, pcedim=5, dst=dst, nua=nua, nub=nub)
            xmcl.append(xmc)
            mccdfl.append(cdfmcx)
            # accytns += cyts
        mckmmed, mckmerrs, mcxc = compmaxdiff(xmcl, mccdfl, xtrth, cdfxtrth)
        print(f'Kolmometer: mc16: {mckmmed:.5f} -- median out of {smpls}')

        xmclt, mccdflt = [], []
        for kkk in range(smpls):
            cyts = np.load(dataprfx + '_pce5_pod8_bfmc32_runs10' + \
                           f'_run{kkk+1}of10_ysoltns.npy')
            xmct, cdfmcxt = comp_cdf(cyts, pcedim=5, dst=dst, nua=nua, nub=nub)
            xmclt.append(xmct)
            mccdflt.append(cdfmcxt)
            # accytns += cyts
        mckmmedt, _, _ = compmaxdiff(xmclt, mccdflt, xtrth, cdfxtrth)
        print(f'Kolmometer: mc32: {mckmmed:.5f} -- median out of {smpls}')

    plt.figure(3)
    # plt.plot(iabsc, np.abs(icdfpcetwo-icdftrth), label='PCE[2]')
    clrs = []
    pltsmpls = int(smpls/2)
    for kkk in range(pltsmpls+1):  # +1 for the legend dummy plot
        clrs.extend([.6])
        clrs.extend([.3])
    clrs.extend([.9])
    # clrs = [.9, .6, .3]
    plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.plasma(clrs))

    for kkk in range(pltsmpls):
        plt.plot(rbxc[kkk][0], rbxc[kkk][1], alpha=.4)
        plt.plot(mcxc[kkk][0], mcxc[kkk][1], alpha=.4)
    plt.plot(np.NaN, np.NaN, label='RB16')
    plt.plot(np.NaN, np.NaN, label='MC16')
    plt.plot(ppxc[0][0], ppxc[0][1], label='pcePOD16')

    plt.title(dst)
    plt.legend()

    plt.show()
