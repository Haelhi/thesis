# Compute quartics --> enumerate_quartics.sage

R = PolynomialRing(QQ,"x")
fields_F = enumerate_totallyreal_fields_prim(4,10^10)
print(fields_F)