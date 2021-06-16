load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

load('functions.sage')
load('data_subfields.sage')

# Kr_list = list_Kr_1
# Kr_list = list_Kr_2
Kr_list = list_Kr_2_not

Kr_list_short = []
for i in range(3,4):
    Kr_list_short.append(Kr_list[i])
Kr_list = Kr_list_short

K_list = CM_one_sextic_from_Kr(Kr_list)

print(K_list)