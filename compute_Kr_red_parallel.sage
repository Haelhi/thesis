load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
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
        if check[0]:
            return((check[0],l))

o = open('output_d8_prime.sage','a')
o.write('K_d8_prime = [')
o.close()

data = []
for i in range(447,len(Kr_d8_prime)):
    data.append(Kr_d8_prime[i])

# data = Kr_d9_prime
chunks = divide_into_chunks(data,len(data))

for x in parallel_reduce(chunks):
    check = list(x)[1]
    print(check)
    if check != None:
        o = open('output_d8_prime.sage','a')
        o.write(str(list(x)[1][1]))
        o.write(',')
        o.close()
    
o = open('output_d8_prime.sage','a')
o.write(']')
o.write('\n\n\n')
o.close()