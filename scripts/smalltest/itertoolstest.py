import numpy as np
import itertools

att = [1, 2, 3]
lenatt = 3*3*3*3
halflen = np.int(np.floor(lenatt/2))

print('first half')
for kkk in itertools.islice(itertools.product(att, repeat=4), 0, halflen):
    print(kkk)

print('second half')
for kkk in itertools.islice(itertools.product(att, repeat=4), halflen, lenatt):
    print(kkk)
