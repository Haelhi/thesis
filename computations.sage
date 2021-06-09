load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations

load('functions.sage')
load('data_subfields.sage')

Kr_list = list_Kr_2

K_list = CM_one_sextic_from_Kr(Kr_list)
print(len(K_list))
print(K_list)