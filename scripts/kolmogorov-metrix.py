import numpy as np
import matplotlib.pyplot as plt
from tikzplotlib import save

import gen_pod_uq.mc_pce_utils as mpu

from kmg_mtrx_helpers import get_evay


def comp_cdf(ysoltens, nua=0., nub=1., dst='beta-2-5', nsample=1000000,
             pcedim=5):
    
    abscissae, _, _, _ = mpu.\
        setup_pce(distribution=dst,
                  distrpars=dict(a=nua, b=nub),
                  pcedim=pcedim, uncdims=4)
    
    evay = get_evay(ysoltens, abscissae)
    
    getsample = mpu.get_nu_sample(distribution=dst,
                                  uncdims=4, nulb=nua, nuub=nub)

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
    nua, nub = 5e-4, 10e-4
    smplsforcdf = int(1e6)

    onlyplots = True
    onlyplots = False

    if onlyplots:
        runs = 1
        smpls = 5  # number of samples for the MC/wRB bases
        pltsmpls = smpls
    else:
        pltsmpls = int(smpls/2)
    pltfilter = 4  # only plot every 4-th data point

    dstl = ['beta-2-5', 'uniform']
    for pltadd, dst in enumerate(dstl):

        Nndstr = f'N12nu{nua:.2e}--{nub:.2e}' + dst
        dataprfx = 'mh-data/cached-data/'  + Nndstr

        ccdfdct = dict(nua=nua, nub=nub, pcedim=5, dst=dst, nsample=smplsforcdf)

        for rrr in range(runs):
            yts = dataprfx + '_pce5_ysoltns.npy'
            ysoltens = np.load(yts)
            xtrth, cdfxtrth = comp_cdf(ysoltens, **ccdfdct)
            # jmin, jmax = xtrth[0], xtrth[-1]

            if onlyplots:
                pass
            else:
                yts = dataprfx + '_pce2_ysoltns.npy'
                ysoltens = np.load(yts)
                xpcetwo, cdfpcetwo = comp_cdf(ysoltens, pcedim=2, dst=dst,
                                              nua=nua, nub=nub, nsample=smplsforcdf)
                ppkmmed, _, ppxc \
                    = compmaxdiff([xpcetwo], [cdfpcetwo], xtrth, cdfxtrth)
                print(f'Kolmometer: {dst}: pce[2]: {ppkmmed:.5f}')

            yts = dataprfx + '_pce5_pod8_bfpce2_run1of1_ysoltns.npy'
            ysoltens = np.load(yts)
            xpodpcef, cdfpodpcef = comp_cdf(ysoltens, **ccdfdct)
            # jmin, jmax = max(jmin, xpodpcef[0]), min(jmax, xpodpcef[-1])
            ppkmmed, _, ppxc \
                = compmaxdiff([xpodpcef], [cdfpodpcef], xtrth, cdfxtrth)
            print(f'Kolmometer: {dst}: pce-16: {ppkmmed:.5f}')

            # accytns = 0
            xrbl, rbcdfl = [], []
            for kkk in range(smpls):
                cyts = np.load(dataprfx + '_pce5_pod8_bfrb_random16_runs10' + \
                               f'_run{kkk+1}of10_ysoltns.npy')
                xrb, cdfrbx = comp_cdf(cyts, **ccdfdct)
                xrbl.append(xrb)
                rbcdfl.append(cdfrbx)
                # accytns += cyts

            rbkmmed, rbkmerrs, rbxc = compmaxdiff(xrbl, rbcdfl, xtrth, cdfxtrth)
            print(f'Kolmometer: {dst}: rb16: {rbkmmed:.5f} -- median out of {smpls}')

            if onlyplots:
                pass
            else:
                xrblt, rbcdflt = [], []
                for kkk in range(smpls):
                    cyts = np.load(dataprfx + '_pce5_pod8_bfrb_random32_runs10' + \
                                   f'_run{kkk+1}of10_ysoltns.npy')
                    xrbt, cdfrbxt = comp_cdf(cyts, **ccdfdct)
                    xrblt.append(xrbt)
                    rbcdflt.append(cdfrbxt)
                    # accytns += cyts

                rbkmmedt, _, _ = compmaxdiff(xrblt, rbcdflt, xtrth, cdfxtrth)
                print(f'Kolmometer: {dst}: rb32: {rbkmmedt:.5f} -- median out of {smpls}')

            xmcl, mccdfl = [], []
            for kkk in range(smpls):
                cyts = np.load(dataprfx + '_pce5_pod8_bfmc16_runs10' + \
                               f'_run{kkk+1}of10_ysoltns.npy')
                xmc, cdfmcx = comp_cdf(cyts, **ccdfdct)
                xmcl.append(xmc)
                mccdfl.append(cdfmcx)

            mckmmed, mckmerrs, mcxc = compmaxdiff(xmcl, mccdfl, xtrth, cdfxtrth)
            print(f'Kolmometer: {dst}: mc16: {mckmmed:.5f} -- median out of {smpls}')

            if onlyplots:
                pass
            else:
                xmclt, mccdflt = [], []
                for kkk in range(smpls):
                    cyts = np.load(dataprfx + '_pce5_pod8_bfmc32_runs10' + \
                                   f'_run{kkk+1}of10_ysoltns.npy')
                    xmct, cdfmcxt = comp_cdf(cyts, **ccdfdct)
                    xmclt.append(xmct)
                    mccdflt.append(cdfmcxt)
                mckmmedt, _, _ = compmaxdiff(xmclt, mccdflt, xtrth, cdfxtrth)
                print(f'Kolmometer: {dst}: mc32: {mckmmed:.5f} -- median out of {smpls}')

        plt.figure(330+pltadd)
        clrs = []
        pltsmpls = int(smpls/2)
        for kkk in range(pltsmpls+1):  # +1 for the legend dummy plot
            clrs.extend([.6])
            clrs.extend([.3])
        clrs.extend([.9])
        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.plasma(clrs))

        for kkk in range(pltsmpls):
            plt.plot(rbxc[kkk][0][::pltfilter], rbxc[kkk][1][::pltfilter], alpha=.4)
            plt.plot(mcxc[kkk][0][::pltfilter], mcxc[kkk][1][::pltfilter], alpha=.4)

        # one dummy point plot to have the labels in full color
        plt.plot(0.65, 0., label='RB16')
        plt.plot(0.65, 0., label='MC16')

        plt.plot(ppxc[0][0][::pltfilter], ppxc[0][1][::pltfilter], label='pcePOD16')

        plt.title(dst)
        plt.legend()
        save('kolmomotor'+dst+'.tikz')

    plt.show()
