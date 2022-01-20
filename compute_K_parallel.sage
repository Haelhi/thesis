pari.allocatemem(80000000000)

load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
load('functions.sage')
load('output_Kr_tk2_q7.sage')
load('output_Kr_tk2_q8.sage')


def CM_one_sextic_from_Kr_print(Kr_list):
    K_list = []
    Kr_pol = Kr_list
    print(Kr_pol[0])
    Kr = CM_Field(Kr_pol[0])
    clKr = Kr.class_group(False)
    gens = clKr.gens()
    rep_gens = []
    for I in gens:
        M = Kr.minkowski_bound().numerical_approx().ceil()
        Ip = I.representative_prime(norm_bound=M)
        rep_gens.append(Ip)
    Phir_set = Kr.CM_types(equivalence_classes=True)
    i = 0
    for Phir in Phir_set:
        i = i + 1
        if test_CM_cl_nr_one_with_class_group(rep_gens, Phir) == True:
            K = Phir.reflex_field()
            print(K.polynomial())
            if K.g() == 3:
                print([pari.polredabs(K.polynomial()),Kr_pol[0],Kr_pol[1],Kr_pol[2]], ',')
                K_list.append([pari.polredabs(K.polynomial()),Kr_pol[0],Kr_pol[1],Kr_pol[2]])
                break
        if i == 6:
            break
    return K_list


# INPUT: list of chunks (output of divide_into_chunks)
# OUTPUT: parallel output with K
@parallel(32)
def parallel_comp(list_of_chunks):
    for chunk in list_of_chunks:
        try:
            CM_one_sextic_from_Kr_print(chunk)
        except Exception:
            print('Error in', chunk)

data = Kr_tk2_q8[100:len(Kr_tk2_18)]
chunks = divide_into_chunks(data,len(data))
list(parallel_comp(chunks))
    
