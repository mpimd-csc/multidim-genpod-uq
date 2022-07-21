import argparse
import json

mesh = 12
nustr = '5.0e-04to1.0e-03'
dist = 'beta-2-5'
dist = 'uniform'
pcelist = [2, 3, 4, 5]

dtfl = 'cached-data/computed-moments.json'

prsr = argparse.ArgumentParser()
prsr.add_argument("--distribution", type=str, help="type of distribution",
                  default=dist)
args = prsr.parse_args()
adist = args.distribution

with open(dtfl, 'r') as fjs:
    dtdct = json.load(fjs)

cddct = dtdct[f'{mesh}'][nustr][adist]

expvlist = [cddct[f'{pcd}']['expv'] for pcd in pcelist]
vrnclist = [cddct[f'{pcd}']['vrnc'] for pcd in pcelist]
ee = [cexp - expvlist[-1] for cexp in expvlist]
ev = [cvrn - vrnclist[-1] for cvrn in vrnclist]

print('**** distribution: ' + adist + '*************')
for kkk, pce in enumerate(pcelist):
    print(f'{pce}: ${expvlist[kkk]:.10f}$/${vrnclist[kkk]:.10f}$ &' +
          f'{ee[kkk]:.1e}/{ev[kkk]:.1e}')
