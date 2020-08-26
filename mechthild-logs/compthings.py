pcedims = [2, 3, 4, 5]
pcevals = [0.88098226, 0.88109208, 0.88101510, 0.88102114]
pceyys = []

for k in range(len(pcevals)-1):
    print('pce{0}-pce{1}: {2}'.format(pcedims[k], pcedims[k+1],
                                      pcevals[k]-pcevals[k+1]))

print('pce{0}-pce{1}: {2}'.format(pcedims[0], pcedims[-1],
                                  pcevals[0]-pcevals[-1]))

for k in range(len(pcevals)):
    print('pce{0}: pceval{1}'.format(pcedims[k], pcevals[k]))
