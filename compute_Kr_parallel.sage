load('Data/data_quartics_A4.sage')
load('Data/k_d5_tk2.sage')
load('Data/k_d5_tk3.sage')
load('functions.sage')
load('elimination_wrt_decomposition.sage')
from datetime import datetime

def prime_check_before_constr(Krplus_pol, dKrplus,dk):
    for prime in list(dKrplus.factor()):
        if dk % prime[0] == 0:
            if elimination_wrt_decomposition_1_3_7_and_15(prime[0], Krplus_pol, dKrplus, dk) == False:
                return False
    return True

def prime_check_after_constr(fKr, fKrplus, dKrplus, dk):
    for prime in list(dKrplus.factor()):
        p = prime[0]
        if dk % p == 0:
            if elimination_wrt_decomposition_19(p, fKrplus, fKr, dKrplus, dk) == False:
                return False
            else:
                pass
        if dk % p != 0:
            if elimination_wrt_decomposition_5_6_and_18(p, fKr, dKrplus, dk) == False:
                return False
            else:
                pass
    return True

def make_list(quadratic,quartic):
    g_list = []
    for k_tup in quadratic:
        for Krplus_tup in quartic:
            g_list.append([k_tup,Krplus_tup])
    return g_list

@parallel(32)
def comp_parallel(list_of_chunks):
    for pair in list_of_chunks:
        k_tup = pair[0]
        print(k_tup)
        Krplus_tup = pair[1]
        construct_Kr_from_k_lmfdb(k_tup,Krplus_tup)

def construct_Kr_from_k_lmfdb(k_tup,Krplus_tup):
    Krplus_pol = Krplus_tup
    Krplus.<a> = NumberField(Krplus_pol)
    dKrplus = Krplus.disc()
    poly_k = k_tup[1]
    dk = k_tup[2]
    k.<a> = NumberField(poly_k)
    hk = k.class_number(False)
    hKrplus = Krplus.class_number(False)
    if prime_check_before_constr(Krplus_pol,dKrplus,dk) == True:
        Kr_rel.<b> = Krplus.extension(poly_k)
        Kr.<c> = Kr_rel.absolute_field()
        hKr = Kr.class_number(False)
        fKr = Kr.polynomial()
        if prime_check_after_constr(fKr, Krplus_pol, dKrplus, dk) == True:
            h = hKr / hKrplus
            lst_h = list((2*h/hk).factor())
            if len(lst_h) == 1:
                if lst_h[0][0] == 2:
                    print([pari.polredabs(Kr.polynomial()),Krplus_pol,poly_k,hKr,hKrplus,hk], ',')

quadratic = k_tk2[965:len(k_tk2)]
quartic = quartic8
prep_list = make_list(quadratic,quartic)
chunks_k = divide_into_chunks(prep_list,len(prep_list))
list(comp_parallel(chunks_k))