# Compute quartics --> enumerate_quartics.sage

L = [1,2,3]
o = open('output_bla.sage','a')
o.write(str(L))
o.close()

R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10^8)
print(fields_F)

o = open('output_quartics_8.sage','a')
o.write(str(fields_F))
o.close()