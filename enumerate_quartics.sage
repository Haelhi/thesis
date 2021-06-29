# Compute quartics --> enumerate_quartics.sage

pari.allocatemem(10000000000)

R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10^9)
print(fields_F)

o = open('output_quartics_9.sage','a')
o.write(str(fields_F))
o.close()