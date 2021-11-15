load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

pari.allocatemem(80000000000)

load('functions.sage')
load('Data/Kr_d7_p_div_dKplus.sage')
load('Data/Kr_d8_p_div_dKplus.sage')
load('Data/Kr_d9_p_div_dKplus.sage')

Kr_list = Kr_d7_prime

n = len(Kr_list)
start = 0
if n < 51:
    end = n
else:
    end = 50
    
while end <= len(Kr_list):
    K_list_partial = CM_one_sextic_from_Kr(Kr_list[start:end])
    o = open('Data/K_CM_prime_d7.sage','a')
    o.write(str(start))
    o.write('\n')
    o.write(str(end))
    o.write('\n')
    o.write(str(K_list_partial))
    o.write('\n\n\n')
    o.close()
    start = start + 50
    end = end + 50
    if end > n:
        K_list_partial = CM_one_sextic_from_Kr(Kr_list[start:n])
        o = open('Data/K_CM_prime_d7.sage','a')
        o.write(str(start))
        o.write('\n')
        o.write(str(end))
        o.write('\n')
        o.write(str(K_list_partial))
        o.write('\n\n\n')
        o.close()