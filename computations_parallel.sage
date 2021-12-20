load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

# pari.allocatemem(80000000000)

load('functions.sage')
load('Data/Kr_d7_p_div_dKplus.sage')
load('Data/Kr_d8_p_div_dKplus.sage')
load('Data/Kr_d9_p_div_dKplus_red.sage')

# INPUT: list of chunks (output of divide_into_chunks)
# OUTPUT: parallel output with K
@parallel(30)
def parallel_comp(list_of_chunks):
    chunks_K = []
    for chunk in list_of_chunks:
        print(chunk)
        K_list_partial = CM_one_sextic_from_Kr_write([chunk])
        chunks_K.append(K_list_partial)
    return chunks_K

print(datetime.now())

o = open('parallel_OUTPUT.sage','a')
o.write('K_d9_prime = [')
o.close()

data = Kr_d9_one
chunks = divide_into_chunks(data,len(data))

for x in parallel_comp(chunks):
    print(list(x))
    
o = open('parallel_output.sage','a')
o.write(']')
o.write('\n\n\n')
o.close()

print(datetime.now())