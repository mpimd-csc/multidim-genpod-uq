import numpy as np

import gen_pod_uq.mc_pce_utils as mpu


def get_evay(ysoltens, abscissae):
    vecy = ysoltens.reshape((1, -1))
    def evay(alphavec):
        psix = np.array(mpu.eva_all_lgrngs(abscissae, alphavec[-1])).\
                reshape((1, -1))
        for calpha in alphavec[:-1]:
            psixn = np.array(mpu.eva_all_lgrngs(abscissae, calpha)).\
                    reshape((1, -1))
            psix = np.kron(psix, psixn)
        return (psix @ vecy.T).item()
    return evay
