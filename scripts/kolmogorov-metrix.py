import numpy as np
import matplotlib.pyplot as plt

import gen_pod_uq.mc_pce_utils as mpu

from kmg_mtrx_helpers import get_evay

def comp_plot_pcf(nua=0., nub=1., nsample=10000):
    nustr = f'nu{nua:.2e}--{nub:.2e}'

    ysltnsstr = 'cached-data/N4' + nustr + dst + '_pce2_ysoltns.npy'
    ysoltens = np.load(ysltnsstr)
    
    abscissae, _, compexpv, _ = mpu.\
        setup_pce(distribution=dst,
                  distrpars=dict(a=nua, b=nub),
                  pcedim=2, uncdims=4)
    
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
    
    plt.figure(1)
    plt.plot(srtdsmpllist, cpfvals)
    plt.show()

if __name__ == '__main__':
    nua, nub = 5e-4, 10e-4
    numean = .5*(nua+nub)
    dst = 'beta-2-5'
    comp_plot_pcf(nua=nua, nub=nub)
