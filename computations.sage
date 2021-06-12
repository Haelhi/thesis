load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

load('functions.sage')
load('data_constructed_Kr.sage')

# Kr_list = list_Kr_1
Kr_list = list_Kr_2

print(len(Kr_list))

Kr_list_short = []
for i in range(0,10):
    Kr_list_short.append(Kr_list[i])
Kr_list = Kr_list_short

K_list = CM_one_sextic_from_Kr(Kr_list)
print(len(K_list))
print(K_list)