from multidim_galerkin_pod.chaos_expansion_utils import get_weighted_gaussqr

absc, weights = get_weighted_gaussqr(N=3, a=3e-4, b=7e-4,
                                     weightfunction='beta',
                                     wfpdict=dict(alpha=2, beta=5))

print('abscisssae:', absc)
print('weights:', weights)

absc, weights = get_weighted_gaussqr(N=5, a=3e-4, b=7e-4,
                                     weightfunction='beta',
                                     wfpdict=dict(alpha=2, beta=5))
print('abscisssae:', absc)
print('weights:', weights)
