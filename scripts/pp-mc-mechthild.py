import numpy as np

# GQ-MC-HD-142890.out:31:mc:1000/1000: estxy[0]=0.8828239386763104
# GQ-MC-HD-142891.out:31:mc:1000/1000: estxy[0]=0.8804753821461133
# GQ-MC-HD-142892.out:31:mc:1000/1000: estxy[0]=0.8765455721067179
# GQ-MC-HD-142893.out:31:mc:1000/1000: estxy[0]=0.8792293283708384
# GQ-MC-HD-142894.out:31:mc:1000/1000: estxy[0]=0.8836003894252067
# GQ-MC-HD-142913.out:31:mc:1000/1000: estxy[0]=0.8809124472623131
# GQ-MC-HD-142914.out:31:mc:1000/1000: estxy[0]=0.8812300016316496
# GQ-MC-HD-142915.out:31:mc:1000/1000: estxy[0]=0.8804289552051776
# GQ-MC-HD-142916.out:31:mc:1000/1000: estxy[0]=0.8850765693327264
# GQ-MC-HD-142917.out:31:mc:1000/1000: estxy[0]=0.8872645563923877


mcreslist = [0.8828239386763104,
             0.8804753821461133,
             0.8765455721067179,
             0.8792293283708384,
             0.8836003894252067,
             0.8809124472623131,
             0.8812300016316496,
             0.8804289552051776,
             0.8850765693327264,
             0.8872645563923877]

mcresarray = np.array(mcreslist)
chunks = np.arange(0, len(mcreslist)+1, 2)
chunklen = 1000
for cks in chunks[1:]:
    nsms = cks*chunklen
    estm = np.mean(mcresarray[:cks])
    print('#sims: {0} -- estimated mean {1:.6f}'.format(nsms, estm))
