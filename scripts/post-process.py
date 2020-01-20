import dolfin_navier_scipy.data_output_utils as dou
import numpy as np
import mat_lib_plots.conv_plot_utils as cpu

jsfstr = 'N3nu3.00e-03--7.00e-03_pcepod_mcpod_bfmc.json'
jsfstr = 'N3nu3.00e-03--7.00e-03_pcepod_mcpod_bfpce.json'
mcref = 1.
pceref = 1.

ddct = dou.load_json_dicts(jsfstr)
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
for timit in tims:
    trntimelist.append(ddct[timit]['traintime'])
    pcepodtimlist = []
    lpcepodreslist = []
    for cpd in poddims:
        pcepodtimlist.append(ddct[timit]['pcepod'][cpd]['elts'])
        cpceres = np.array(ddct[timit]['pcepod'][cpd]['pceres'])
        lpcepodreslist.append(cpceres.flatten())
    teltlist.append(pcepodtimlist)
    pcepodreslist.append(lpcepodreslist)
    crmlist.append(ddct[timit]['comp-redmod-elts'])
    rmprjelist.append(np.array(ddct[timit]['redmod-prj-errs']).flatten())

pcepodresarray = np.array(pcepodreslist)
pceerrarray = pceref - pcepodresarray
print('***pce errrors***')
cpu.print_nparray_tex(np.median(pceerrarray, axis=0),
                      formatit='math', fstr='.4e')

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


# ## MC POD
print('*** MC POD ***')
print('*** with basis from {0} ***'.format(basisfrom))

poddims = list((ddct['0']['mcpod']).keys())
mcruns = ddct['0']['mcpod'][poddims[0]]['mcruns']

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
print('***mc errrors***')
cpu.print_nparray_tex(np.median(mcerrarray, axis=0),
                      formatit='math', fstr='.4e')

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

print('*** mc elts ***')
cpu.print_nparray_tex(teltarray.min(axis=0), formatit='texttt', fstr='.2f')
