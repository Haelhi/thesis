R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10^9)
print(len(fields_F))
list_F = []
for F in fields_F:
    f = R(list(pari.Col(F[1])))
    K.<a> = NumberField(f)
    L.<b> = K.galois_closure()
    if L.degree() == 12:
        print(L.degree()); print(f)
        list_F.append(pari(f).polredabs())
print(list_F)