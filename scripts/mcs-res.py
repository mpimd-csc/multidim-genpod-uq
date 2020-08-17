mcl = [[-5.592134638754766e-05, 1e6],
       [-6.441295571105171e-05, 1e6],
       [-2.7292256647813406e-05, 1e6],
       [-5.738494934280777e-05, 1e6],
       [-3.920474118683737e-05, 1e6]]


accval, accruns = 0, 0
for mcrun in mcl:
    accval += mcrun[0]*mcrun[1]
    accruns += mcrun[1]
    print(accruns, ': ', accval/accruns)
