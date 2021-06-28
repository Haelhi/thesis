# Compute quartics --> enumerate_quartics.sage

R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10000)
# print(fields_F)
o = open('output_quartics.txt','a')
o.write(str(fields_F))
o.close()

"""
R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10000)
#print(len(fields_F))
for F in fields_F:
    f = R(list(pari.Col(F[1])))
    K.<a> = NumberField(f)
    L.<b> = K.galois_closure()
    if L.degree() == 12:
        f = pari(f).polredabs()
        o = open('OUTPUT_quartic_A4.sage','a')  
        o.write(str(f))
        o.write(',\n\n') 
        o.close()
    if L.degree() == 24:
        f = pari(f).polredabs()
        o = open('OUTPUT_quartic_S4.sage','a')  
        o.write(str(f))
        o.write(',\n\n') 
        o.close()
"""