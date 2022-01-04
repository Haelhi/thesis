load('Data/data_quartics_A4.sage')
load('Data/k_d5_tk1.sage')
load('functions.sage')

pari.allocatemem(83613065216)

def check_2hKrstar_div_hk_poweroftwo(poly_Kr, poly_k, hKr, hKrplus, hk):
    hKrstar = hKr / hKrplus
    h = 2*hKrstar / hk
    for p in list(h.factor()):
        if p[0] % 2 == 1:
            return (False,0)
    return (True,h.factor())

def construct_Kr_from_k_lmfdb(f_k,f_Krplus,hk):
    Krplus.<a> = NumberField(f_Krplus)
    hKrplus = Krplus.class_number(False)
    Kr_rel.<b> = Krplus.extension(f_k)
    Kr.<c> = Kr_rel.absolute_field()
    hKr = Kr.class_number(False)
    f_Kr = Kr.polynomial()
    check = check_2hKrstar_div_hk_poweroftwo(f_Kr, f_k, hKr, hKrplus, hk)
    if check[0]:
        return([f_Kr,f_Krplus,f_k,hKr,hKrplus,hk,check[1]])

def make_chunks(k_list,Krplus_list):
    comp_list = []
    for k in k_list:
        for Krplus in Krplus_list:     
            comp_list.append([Krplus,k[1],k[3]])
    return(comp_list)

@parallel(30)
def Kr_from_k_parallel(list_of_chunks):
    chunks_Kr = []
    for poly in list_of_chunks:
        Krplus_poly = poly[0]
        k_poly = poly[1]
        hk = poly[2]
        Kr_values = construct_Kr_from_k_lmfdb(k_poly,Krplus_poly,hk)
        chunks_Kr.append(Kr_values)
    return(chunks_Kr)


o = open('Kr_d7_tk1.sage','a')
o.write('Kr_d7_tk1 = [')
o.close()

quadratic = k_tk1[1000:1500]
quartic = quartic7
list_Krplus_k = make_chunks(quadratic,quartic)
n = len(list_Krplus_k)
chunks = divide_into_chunks(list_Krplus_k,n)

for x in Kr_from_k_parallel(chunks):
    if x[1][0] != None:
        o = open('Kr_d7_tk1.sage','a')
        o.write(str(x[1][0]))
        o.write(',')
        o.close()
    
o = open('Kr_d7_tk1.sage','a')
o.write(']')
o.write('\n\n\n')
o.close()