# Compute quartics --> enumerate_quartics.sage

R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10^9)
print(len(fields_F))
list_A4 = []
list_S4 = []
for F in fields_F:
    f = R(list(pari.Col(F[1])))
    K.<a> = NumberField(f)
    L.<b> = K.galois_closure()
    if L.degree() == 12:
        f = pari(f).polredabs()
        list_A4.append([f,F[0]])
    if L.degree() == 24:
        f = pari(f).polredabs()
        list_S4.append([f,F[0]])
print(list_A4)
print(list_S4)