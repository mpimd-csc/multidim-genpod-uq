from scipy.stats import beta
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

a, b = 2, 5

nulb, nuub = 3e-4, 7e-4

x = np.linspace(nulb, nuub, 100)

ax.plot(x, beta.pdf(x, a, b, loc=nulb, scale=nuub-nulb),
        'r-', lw=5, alpha=0.6, label='beta pdf')
plt.show()
