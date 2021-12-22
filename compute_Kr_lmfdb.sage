load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
load('Data/data_quartics_A4.sage')
load('Data/k_d5_tk1.sage')
load('functions.sage')

def check_2hKrstar_div_hk_poweroftwo(poly_Kr, poly_k, hKr, hKrplus, hk):
    hKrstar = hKr / hKrplus
    h = 2*hKrstar / hk
    for p in list(h.factor()):
        if p[0] % 2 == 1:
            return (False,0)
    return (True,h.factor())

def construct_Kr_from_k_lmfdb(k_list,f_Krplus):
    Krplus.<a> = NumberField(f_Krplus)
    hKrplus = Krplus.class_number(False)
    Kr_list = []
    for i in k_list:
        f_k = i[1]
        hk = i[2]
        Kr_rel.<b> = Krplus.extension(f_k)
        Kr.<c> = Kr_rel.absolute_field()
        hKr = Kr.class_number(False)
        f_Kr = Kr.polynomial()
        check = check_2hKrstar_div_hk_poweroftwo(f_Kr, f_k, hKr, hKrplus, hk)
        if check[0]:
            Kr_list.append([f_Kr,f_Krplus,f_k,hKr,hKrplus,hk,check[1]])
    return(Kr_list)

@parallel(30)
def Kr_from_k_parallel(k_list,list_of_chunks):
    chunks_Kr = []
    for chunk in list_of_chunks:
        print(chunk)
        Kr_list = construct_Kr_from_k_lmfdb(k_list,chunk[0])
        chunks_Kr.append(Kr_list)
    return(chunks_Kr)



o = open('parallel_OUTPUT.sage','a')
o.write('Kr_d7_tk1 = [')
o.close()

quadratic = k_tk1
quartic = quartic7
chunks = divide_into_chunks(quartic,len(quartic))

for x in Kr_from_k_parallel(quadratic,chunks):
    print(list(x))
    
o = open('parallel_OUTPUT.sage','a')
o.write(']')
o.write('\n\n\n')
o.close()