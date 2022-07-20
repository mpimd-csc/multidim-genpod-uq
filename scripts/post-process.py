import argparse
import json
import matplotlib.pyplot as plt
from tikzplotlib import save
import numpy as np

dst = 'uniform'
dst = 'beta-2-5'
nustr = 'nu5.00e-04--1.00e-03'

prsr = argparse.ArgumentParser()
prsr.add_argument("--distribution", type=str, help="type of distribution",
                  default=dst)
args = prsr.parse_args()
adst = args.distribution

pltdict = {'rb': dict(label='wRB', mrkr='.', size=(7, 5), colrs=[.5, .6]),
           'pce': dict(label='PCE', mrkr='.', size=(7, 5), colrs=[.8, .9]),
           'mc': dict(label='Random', mrkr='.', size=(7, 5), colrs=[.2, .3]),
           }

fstrl = [f'cached-data/N12{nustr}{dst}_pcepod_bfmc16_runs10.json',
         f'cached-data/N12{nustr}{dst}_pcepod_bfmc32_runs10.json',
         f'cached-data/N12{nustr}{dst}_pcepod_bfmc64_runs10.json',
         f'cached-data/N12{nustr}{dst}_pcepod_bfrb_random32_runs10.json',
         f'cached-data/N12{nustr}{dst}_pcepod_bfrb_random16_runs10.json',
         f'cached-data/N12{nustr}{dst}_pcepod_bfpce2.json',
         ]

lbllst = ['MC-16', 'MC-32', 'MC-64', 'RB-16', 'RB-32', 'PCE-16']
sumclrs = [.2, .26, .32, .5, .59, .8]

print('distribution: ' + adst)

fignum = 101

for kkk, jsfstr in enumerate(fstrl):

    with open(jsfstr, 'r') as f:
        ddct = json.load(f)

    basisfrom = ddct['0']['basisfrom']
    truthexpy = ddct['truthexpy']
    truthvrnc = ddct['truthvrnc']

    # ## Collect the data

    # ## PCE POD
    print('*** PCE POD ***')
    print('*** with basis from {0} ***'.format(basisfrom))

    poddims = list((ddct['0']['pcepod']).keys())
    pcedims = ddct['0']['pcepod'][poddims[0]]['pcedims']

    tims = list(ddct.keys())

    teltlist = []
    trntimelist = []
    crmlist = []
    rmprjelist = []
    pcepodreslist = []
    pcepodeyyslist = []
    ntims = 0  # counting timings
    for timit in tims:
        try:
            trntimelist.append(ddct[timit]['traintime'])
            pcepodtimlist = []
            lpcepodreslist = []
            lpcepodeyyslist = []
            for cpd in poddims:
                pcepodtimlist.append(ddct[timit]['pcepod'][cpd]['elts'])
                cpceres = np.array(ddct[timit]['pcepod'][cpd]['pceres'])
                cpceyys = np.array(ddct[timit]['pcepod'][cpd]['pcepodeyys'])
                lpcepodreslist.append(cpceres.flatten())
                lpcepodeyyslist.append(cpceyys.flatten())
            teltlist.append(pcepodtimlist)
            pcepodreslist.append(lpcepodreslist)
            pcepodeyyslist.append(lpcepodeyyslist)
            crmlist.append(ddct[timit]['comp-redmod-elts'])
            rmprjelist.append(np.array(ddct[timit]
                                       ['redmod-prj-errs']).flatten())
            ntims += 1
        except TypeError:
            pass  # no timit key

    pcepodresarray = np.array(pcepodreslist)
    pceeyys = np.array(pcepodeyyslist)

    pceerrarray = np.abs(truthexpy - pcepodresarray)/truthexpy
    pcevrnce = np.abs(truthvrnc - (pceeyys - pcepodresarray**2))/truthvrnc

    poddimsints = [np.int(pd) for pd in poddims]
    cpltd = pltdict[basisfrom]

    plt.figure(fignum, figsize=(6, 3))
    plt.rcParams["axes.prop_cycle"] = \
        plt.cycler("color", plt.cm.plasma(cpltd['colrs']))
    if ntims > 1:
        for ct in range(ntims):
            plt.semilogy(poddimsints, pcevrnce[ct, :, -1], cpltd['mrkr'],
                         markersize=cpltd['size'][1], alpha=.1)
            plt.semilogy(poddimsints, pceerrarray[ct, :, -1], cpltd['mrkr'],
                         markersize=cpltd['size'][0], alpha=.1)
    else:
        pass
    plt.semilogy(poddimsints, np.median(pcevrnce[:, :, -1], axis=0),
                 cpltd['mrkr'], markersize=cpltd['size'][1],
                 label=lbllst[kkk]+'--E')
    plt.semilogy(poddimsints, np.median(pceerrarray[:, :, -1], axis=0),
                 cpltd['mrkr'],
                 markersize=cpltd['size'][0],
                 label=lbllst[kkk]+'--E')
    plt.xlabel('ROM dimension')
    plt.ylim([.5*1e-6, 1])
    plt.title('ROM Approximation Errors')
    plt.legend()
    plt.tight_layout()
    save('pics/' + dst+'-'+lbllst[kkk]+'.tex')
    fignum += 1

    plt.figure(10101, figsize=(6, 4))
    plt.rcParams["axes.prop_cycle"] = \
        plt.cycler("color", plt.cm.plasma(sumclrs))
    plt.semilogy(poddimsints, np.median(pceerrarray[:, :, -1], axis=0),
                 label=lbllst[kkk], linewidth=4)

plt.xlabel('ROM dimension')
plt.ylim([.5*1e-6, 1e-1])
plt.title('Median ROM Approximation Errors for the Expected Value')
plt.legend()
plt.tight_layout()
save('pics/'+dst+'-eerrors.tex')

plt.show()
