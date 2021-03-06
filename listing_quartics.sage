import json

with open('cubics_0-6_C3.json') as u:
  listC3 = json.load(u)
with open('cubics_0-6_S3.json') as w:
  listS3 = json.load(w)

Krplus_C3 = []
R = PolynomialRing(QQ,'x')
for i in listC3:    
    j = i[0]
    a3 = ZZ(-j[0]); a2 = ZZ(j[1]); a1 = ZZ(-j[2])
    if a3.is_square() == True:
        k.<b> = NumberField(R(i[0]))
        coeff = [a1^2-4*a2, -8*sqrt(a3), -2*a1, 0, 1]
        f = list(pari.polredabs(R(coeff)))
        K.<a> = NumberField(R(f))
        sK = K.signature()
        if f not in Krplus_C3 and sK == (4,0):
            Krplus_C3.append(f)
print(Krplus_C3)

Krplus_S3 = []
R = PolynomialRing(QQ,'x')
for i in listS3:    
    j = i[0]
    a3 = ZZ(-j[0]); a2 = ZZ(j[1]); a1 = ZZ(-j[2])
    if a3.is_square() == True:
        coeff = [a1^2-4*a2, -8*sqrt(a3), -2*a1, 0, 1]
        f = list(pari.polredabs(R(coeff)))
        K.<a> = NumberField(R(f))
        sK = K.signature()
        if f not in Krplus_S3 and sK == (4,0):
            Krplus_S3.append(f)
print(Krplus_S3)