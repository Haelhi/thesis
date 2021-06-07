#from recip import *
from itertools import combinations

quartic = [
[1, x^4 - x^3 - 7*x^2 + 2*x + 9, 26569, 4, []]]

def construct_Kr(quartic):
    Kr_list = []
    for quart in quartic:
        Krplus.<f0> = NumberField(quart[1])
        R0.<t0> = PolynomialRing(Krplus)
        dKrplus = Krplus.discriminant()
        kprimes = list(dKrplus.factor())
        for i in range(1,len(kprimes)+1):
            check_isom = 0
            if i == 1:
                prime_combos = kprimes
                for p in prime_combos:
                    k.<x> = QuadraticField(-p[0])
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
                        Kr_list.append([Kr.polynomial(),k.polynomial(),Krplus.class_number()])
            else:
                prime_combos = list(combinations(kprimes,i))
                for p in prime_combos:
                    d = 1
                    for n in range(1,i+1):
                        d = d * p[i-1][0]
                    k.<x> = QuadraticField(-d)
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
                        Kr_list.append([Kr.polynomial(),k.polynomial(),Krplus.class_number()])
    return Kr_list

Kr_list = construct_Kr(quartic)
print(Kr_list)