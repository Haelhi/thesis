from recip import *
from datetime import datetime 

load('functions.sage')
load('Data/Kr_d7_p_div_dKplus.sage')
load('Data/Kr_d8_p_div_dKplus.sage')
load('Data/Kr_d9_p_div_dKplus.sage')

# INPUT: list of chunks (output of divide_into_chunks)
# OUTPUT: parallel output with Kr reduced
@parallel(30)
def parallel_reduce(list_of_chunks):
    list_Kr = []
    for l in list_of_chunks:
        Kr.<a> = NumberField(l[0])
        check = check_indices_poweroftwo(Kr,l[2])
        print(check)
        if check == True:
            print(l)
            list_Kr.append(l)
    return list_Kr

o = open('output_Kr_red.sage','a')
o.write('K_d8_prime = [')
o.close()

data = Kr_d8_prime
chunks = divide_into_chunks(data,len(Kr_d8_prime))

for x in parallel_reduce(chunks):
    print(list(x))
    
o = open('output_Kr_red.sage','a')
o.write(']')
o.write('\n\n\n')
o.close()