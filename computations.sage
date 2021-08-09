load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

load('functions.sage')
# load('quartic8.sage')
load('data_subfields.sage')

#quartic = quartic8
#Kr_list0 = construct_Kr(quartic)
#Kr_list = []
#for Kr in Kr_list0:
#    if Kr not in Kr_list:
#        Kr_list.append(Kr)
#print(Kr_list)

Kr_list = list_Kr_d8

Kr_list_short = []
for i in range(200,250):
    Kr_list_short.append(Kr_list[i])
Kr_list = Kr_list_short

K_list = CM_one_sextic_from_Kr(Kr_list)

print(K_list)