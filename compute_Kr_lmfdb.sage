load('Data/data_quartics_A4.sage')
load('Data/k_d5_tk1.sage')
load('functions.sage')

pari.allocatemem(83613065216)

def construct_Kr_from_k_lmfdb(f_k,f_Krplus,hk):
    Krplus.<a> = NumberField(f_Krplus)
    hKrplus = Krplus.class_number(False)
    Kr_rel.<b> = Krplus.extension(f_k)
    Kr.<c> = Kr_rel.absolute_field()
    hKr = Kr.class_number(False)
    f_Kr = Kr.polynomial()
    check = check_2hKrstar_div_hk_poweroftwo(f_Kr, f_k, hKr, hKrplus, hk)
    if check[0]:
        return([pari.polredabs(Kr.polynomial()),f_Krplus,f_k,hKr,hKrplus,hk,check[1]])

quadratic = k_tk1[1500:2000]
quartic = quartic7

o = open('output_Kr_lmfdb.sage', 'a')
o.write('[')
o.close()

for j in quadratic:
    for i in quartic:
        Kr_data = construct_Kr_from_k_lmfdb(j[1],i,j[3])
        if Kr_data != None:
            o = open('output_Kr_lmfdb.sage', 'a')
            o.write(str(Kr_data))
            o.write(',')
            o.close()

o = open('output_Kr_lmfdb.sage', 'a')
o.write('0]')
o.close()