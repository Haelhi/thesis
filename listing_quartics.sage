load("cubics_6-7.sage")

def check_hessian(hess):
    for i in hess:
        if hess[0] != i and hess[0] != -i:
            return False
    return True
    
Kplus_C3 = [i for i in cubics if check_hessian(i[1])==True]
Kplus_S3 = [i for i in cubics if check_hessian(i[1])==False]

o = open('./output/cubics_6-7_C3.sage','a')
o.write(str(Kplus_C3))
o.close()
o = open('./output/cubics_6-7_S3.sage','a')
o.write(str(Kplus_S3))
o.close()