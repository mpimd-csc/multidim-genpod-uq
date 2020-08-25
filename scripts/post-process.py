import dolfin_navier_scipy.data_output_utils as dou
import numpy as np
import mat_lib_plots.conv_plot_utils as cpu

N = 3
podpcebas = 2
filedir = ''
# filedir = '../mechthild-logs/mechthildlogs/'
jsfstr = 'N{0}nu3.00e-04--7.00e-04_pcepod{1}_bfpce.json'.format(N, podpcebas)
# jsfstr = 'N{0}nu3.00e-04--7.00e-04_pcepod_mcpod_bfpce.json'.format(N)
# jsfstr = 'N{0}nu3.00e-04--7.00e-04_pcepod_mcpod_bfmc.json'.format(N)
mcref = 0.881759
pceref = 0.88159
pceeyyref = 0.
mcref = pceref

mcpod = False

ddct = dou.load_json_dicts(filedir+jsfstr)
basisfrom = ddct['0']['basisfrom']

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
for timit in tims:
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
    rmprjelist.append(np.array(ddct[timit]['redmod-prj-errs']).flatten())

pcepodresarray = np.array(pcepodreslist)
pceerrarray = pceref - pcepodresarray

pcepodeyyarray = np.array(pcepodeyyslist)
pcepodvrncserrarray = pcepodeyyarray - pcepodeyyarray**2 -\
        (pceeyyref - pceref**2)

if basisfrom == 'pce':
    trainpceexpv = ddct['0']['training-pce-expv']
    print('***training errror***')
    print('{0:.4e}'.format(pceref-trainpceexpv[0]))

print(poddims)
print(pcedims)
print('***pce errrors E(y)***')
cpu.print_nparray_tex(np.median(pceerrarray, axis=0),
                      formatit='math', fstr='.2e')

print('***pce errrors V(y)***')
cpu.print_nparray_tex(np.median(pcepodvrncserrarray, axis=0),
                      formatit='math', fstr='.2e')

print('*** training time (min out of {0})***'.format(len(tims)))
cpu.print_nparray_tex((np.array(trntimelist)).min(),
                      formatit='texttt', fstr='.2f')

print('*** poddims and ' +
      'comp red mod (min out of {0})***'.format(len(tims)))
print(poddims)
cpu.print_nparray_tex((np.array(crmlist)).min(axis=0),
                      formatit='texttt', fstr='.2f')
print('*** comp red projection error (med out of {0})***'.format(len(tims)))
print(np.median(np.array(rmprjelist), axis=0))

teltarray = np.array(teltlist)

print('*** pce elts ***')
cpu.print_nparray_tex(teltarray.min(axis=0), formatit='texttt', fstr='.2f')


if mcpod:
    # ## MC POD
    print('*** MC POD ***')
    print('*** with basis from {0} ***'.format(basisfrom))

    poddims = list((ddct['0']['mcpod']).keys())
    mcruns = ddct['0']['mcpod'][poddims[0]]['mcruns']
    print('*** mc runs ***')
    print(mcruns)

    tims = list(ddct.keys())

    teltlist = []
    trntimelist = []
    crmlist = []
    rmprjelist = []
    mcpodreslist = []
    for timit in tims:
        trntimelist.append(ddct[timit]['traintime'])
        mcpodtimlist = []
        lmcpodreslist = []
        for cpd in poddims:
            mcpodtimlist.append(ddct[timit]['mcpod'][cpd]['elt'])
            cmcres = np.array(ddct[timit]['mcpod'][cpd]['mcres'])
            lmcpodreslist.append(cmcres.flatten())
        teltlist.append(mcpodtimlist)
        mcpodreslist.append(lmcpodreslist)
        crmlist.append(ddct[timit]['comp-redmod-elts'])
        rmprjelist.append(np.array(ddct[timit]['redmod-prj-errs']).flatten())

    mcpodresarray = np.array(mcpodreslist)
    mcerrarray = mcref - mcpodresarray
    print('***mc errrors median (out of {0})***'.format(len(tims)))
    cpu.print_nparray_tex(np.median(mcerrarray, axis=0),
                          formatit='math', fstr='.4e')

    print('***mc errrors min (out of {0})***'.format(len(tims)))
    cpu.print_nparray_tex(np.min(mcerrarray, axis=0),
                          formatit='math', fstr='.4e')

    print('*** training time (min out of {0})***'.format(len(tims)))
    cpu.print_nparray_tex((np.array(trntimelist)).min(),
                          formatit='texttt', fstr='.2f')

    print('*** poddims and ' +
          'comp red mod (min out of {0})***'.format(len(tims)))
    print(poddims)
    cpu.print_nparray_tex((np.array(crmlist)).min(axis=0),
                          formatit='texttt', fstr='.2f')
    print('*** comp red projection error (med out of {0})***'.
          format(len(tims)))
    print(np.median(np.array(rmprjelist), axis=0))

    teltarray = np.array(teltlist)

    print('*** mc elts ***')
    cpu.print_nparray_tex(teltarray.min(axis=0), formatit='texttt', fstr='.2f')
