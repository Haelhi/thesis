load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
load('Data/Kr_d8_p_div_dKplus.sage')

from itertools import combinations
from datetime import datetime

pari.allocatemem(83613065216)

load('functions.sage')

Kr_list_one = Kr_d8_one
Kr_list_prime = Kr_d8_prime

Kr_list_red_one = []
Kr_list_red_prime = []

for l in Kr_list_one:
    Kr.<a> = NumberField(l[0])
    check = check_indices_poweroftwo(Kr)
    print(check)
    if check == True:
        Kr_list_red_one.append(l)
        
for l in Kr_list_prime:
    Kr.<a> = NumberField(l[0])
    check = check_indices_poweroftwo(Kr)
    print(check)
    if check == True:
        Kr_list_red_prime.append(l)
        
print(Kr_list_red_one)
print('---------------')
print(Kr_list_red_prime)