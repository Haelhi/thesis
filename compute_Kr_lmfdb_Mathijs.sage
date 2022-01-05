load('Data/data_quartics_A4.sage')
load('Data/k_d5_tk1.sage')
load('functions.sage')
load('elimination_wrt_decomposition.sage')

quadratic = k_tk1
quartic = quartic7

def prime_check_before_constr(dKrplus,dk):
    for prime in list(dKrplus.factor()):
        if dk % prime[0] == 0:
            if elimination_wrt_decomposition_1_3_7_and_15(prime[0], Krplus_pol, dKrplus, dk) == False:
                return False
    return True

parts_list = []

for k_tup in quadratic[500:1000]:
    print(k_tup[1])
    for Krplus_pol in quartic:
        Krplus.<a> = NumberField(Krplus_pol)
        dKrplus = Krplus.disc()
        poly_k = k_tup[1]
        dk = k_tup[2]
        hk = k_tup[3]
        hKrplus = Krplus.class_number(False)
        if prime_check_before_constr(dKrplus,dk) == True:
            Kr_rel.<b> = Krplus.extension(poly_k)
            Kr.<c> = Kr_rel.absolute_field()
            hKr = Kr.class_number(False)
            f_Kr = Kr.polynomial()
            h = hKr / hKrplus
            lst_h = list((2*h/hk).factor())
            if len(lst_h) == 1:
                if lst_h[0][0] == 2:
                    print('Yes')
                    parts_list.append([pari.polredabs(Kr.polynomial()),Krplus_pol,poly_k,hKr,hKrplus,hk])
print(parts_list)