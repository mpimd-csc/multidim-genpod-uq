import numpy as np

pcelist = [2, 3, 4, 5]
pces = [0.8809823,
        0.8810921,
        0.8810151,
        0.8810211]

for pcers in pces:
    print(pcers-pces[-1])

# GQ-MC-HD-142890.out:mc:1000/1000: estxy[0]=0.8828239386763104
# GQ-MC-HD-142891.out:mc:1000/1000: estxy[0]=0.8804753821461133
# GQ-MC-HD-142892.out:mc:1000/1000: estxy[0]=0.8765455721067179
# GQ-MC-HD-142893.out:mc:1000/1000: estxy[0]=0.8792293283708384
# GQ-MC-HD-142894.out:mc:1000/1000: estxy[0]=0.8836003894252067
# GQ-MC-HD-142913.out:mc:1000/1000: estxy[0]=0.8809124472623131
# GQ-MC-HD-142914.out:mc:1000/1000: estxy[0]=0.8812300016316496
# GQ-MC-HD-142915.out:mc:1000/1000: estxy[0]=0.8804289552051776
# GQ-MC-HD-142916.out:mc:1000/1000: estxy[0]=0.8850765693327264
# GQ-MC-HD-142917.out:mc:1000/1000: estxy[0]=0.8872645563923877
# GQ-MC-HD-143002.out:mc:1000/1000: estxy[0]=0.8803112620439154
# GQ-MC-HD-143003.out:mc:1000/1000: estxy[0]=0.875519806999519
# GQ-MC-HD-143004.out:mc:1000/1000: estxy[0]=0.8881222875035939
# GQ-MC-HD-143005.out:mc:1000/1000: estxy[0]=0.8817222086354587
# GQ-MC-HD-143006.out:mc:1000/1000: estxy[0]=0.8806005056808069

mcreslist = [0.8828239386763104,
             0.8804753821461133,
             0.8765455721067179,
             0.8792293283708384,
             0.8836003894252067,
             0.8809124472623131,
             0.8812300016316496,
             0.8804289552051776,
             0.8850765693327264,
             0.8872645563923877,
             0.8803112620439154,
             0.875519806999519,
             0.8881222875035939,
             0.8817222086354587,
             0.8806005056808069]

chunkl = 5
npit = 1000
mcrsar = np.array(mcreslist)
chunks = np.arange(0, len(mcreslist), chunkl)
print(chunks)

for chk in chunks[1:]:
    estk = np.mean(mcrsar[:chk])
    print('{0} mc runs: estimated mean {1}'.format((chk)*npit, estk))

estk = np.mean(mcrsar)
print('{0} mc runs: estimated mean {1}'.format(len(mcreslist)*npit, estk))
