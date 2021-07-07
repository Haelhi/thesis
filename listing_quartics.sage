load('cubics_6-7_C3.sage')
load('cubics_6-7_S3.sage')

def check_hessian(hess):
    for i in hess:
        if hess[0] != i and hess[0] != -i:
            return False
    return True
    
listC3 = [i for i in list_entries if check_hessian(i[1])==True]
listS3 = [i for i in list_entries if check_hessian(i[1])==False]

Krplus_C3 = []
R = PolynomialRing(QQ,'x')
for i in listC3:    
    j = list(pari.polredabs(R(i[0])))
    a3 = ZZ(-j[0]); a2 = ZZ(j[1]); a1 = ZZ(-j[2])
    if a3.is_square() == True:
        k.<b> = NumberField(R(i[0]))
        coeff = []
        coeff.append(a1^2-4*a2)
        coeff.append(-8*sqrt(a3))
        coeff.append(-2*a1)
        coeff.append(0)
        coeff.append(1)
        f = list(pari.polredabs(R(coeff)))
        K.<a> = NumberField(R(f))
        sK = K.signature()
        if f not in Krplus_C3 and sK == (4,0):
            Krplus_C3.append(f)
print(Krplus_C3)

Krplus_S3 = []
R = PolynomialRing(QQ,'x')
for i in listS3:    
    j = list(pari.polredabs(R(i[0])))
    a3 = ZZ(-j[0]); a2 = ZZ(j[1]); a1 = ZZ(-j[2])
    if a3.is_square() == True:
        coeff = []
        coeff.append(a1^2-4*a2)
        coeff.append(-8*sqrt(a3))
        coeff.append(-2*a1)
        coeff.append(0)
        coeff.append(1)
        f = list(pari.polredabs(R(coeff)))
        K.<a> = NumberField(R(f))
        sK = K.signature()
        if f not in Krplus_S3 and sK == (4,0):
            Krplus_S3.append(f)
print(Krplus_S3)