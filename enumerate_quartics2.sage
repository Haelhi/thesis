# Compute quartics --> enumerate_quartics.sage

fields_F = enumerate_totallyreal_fields_prim(4,10^8)
print(fields_F)

o = open('output_quartics_8.sage','a')
o.write(str(fields_F))
o.close()