from recip import *

K.<a> = NumberField(x^6 + 10*x^4 + 21*x^2 + 11)
Kplus.<g> = NumberField(x^3 - 10*x^2 + 21*x - 11)
L.<b> = K.galois_closure()
Lplus.<c> = NumberField(x^12 + 168*x^11 + 11536*x^10 + 422400*x^9 + 9090208*x^8 + 120397568*x^7 + 1003543168*x^6 + 5266527744*x^5 + 16963399936*x^4 + 31521900544*x^3 + 29859194880*x^2 + 11002355712*x + 2985984)
Kr.<d> = NumberField(x^8 + 62*x^6 + 881*x^4 + 3339*x^2 + 1600)
k.<e> = NumberField(x^2 + 16*x + 5883)
Krplus.<f> = NumberField(x^4 - 62*x^3 + 881*x^2 - 3339*x + 1600)
A.<g> = NumberField(x^6 + 9900*x^4 + 19131552*x^2 + 7513995456)
# nfs = []
# for i in L.subfields():
#     f = i[0].absolute_polynomial()
#     if f.degree() == 12:
#         nfs.append(i[0])
# for j in nfs:
#     if len(j.real_embeddings()) == 12:
#         print(j)
dk = k.discriminant().factor()
dLplus = Lplus.discriminant().factor()
dL = L.discriminant().factor()
dK = K.discriminant().factor()
dKplus = Kplus.discriminant().factor()
dKr = Kr.discriminant().factor()
dKrplus = Krplus.discriminant().factor()
print(dk)
print(dK/dKplus^2)
# print(dKrplus)
# print(dKr/dk^2)
A.ideal(11).factor()
G = L.galois_group()
H = G.decomposition_group(L.primes_above(11)[0])
H.fixed_field()[0].is_isomorphic(K)
# E = G.inertia_group(L.primes_above(37)[0])
# E.fixed_field()









