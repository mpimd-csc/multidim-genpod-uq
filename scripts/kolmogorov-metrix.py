import numpy as np

import gen_pod_uq.mc_pce_utils as mpu


nua, nub = 5e-4, 10e-4
numean = .5*(nua+nub)
dst = 'beta-2-5'
nustr = f'nu{nua:.2e}--{nub:.2e}'

ysltnsstr = 'cached-data/N4' + nustr + dst + '_pce2_ysoltns.npy'
ysoltens = np.load(ysltnsstr)

abscissae, _, compexpv, _ = mpu.\
    setup_pce(distribution=dst,
              distrpars=dict(a=nua, b=nub),
              pcedim=2, uncdims=4)

expv = compexpv(ysoltens)
vecy = ysoltens.reshape((1, -1))

print(f'expv: {expv.item():.4e}')

def evay(alphavec):
    psix = np.array(mpu.eva_all_lgrngs(abscissae, alphavec[-1])).\
            reshape((1, -1))
    for calpha in alphavec[:-1]:
        psixn = np.array(mpu.eva_all_lgrngs(abscissae, calpha)).\
                reshape((1, -1))
        psix = np.kron(psix, psixn)
    return (psix @ vecy.T).item()

yamn = evay([numean]*4)
print(f'y(amean): {yamn:.4e}')
yamn = evay([nua]*4)
print(f'y(amin): {yamn:.4e}')
yamn = evay([nub]*4)
print(f'y(amax): {yamn:.4e}')

getsample = mpu.get_nu_sample(distribution=dst,
                              uncdims=4, nulb=nua, nuub=nub)
randa = getsample(1)
yamn = evay(randa.flatten())
print(f'y(arnd): {yamn:.4e}')
