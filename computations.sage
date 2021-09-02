load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

pari.allocatemem(83613065216)

load('functions.sage')
load('data_quartics_S4.sage')
load('data_Kr_d10.sage')

#quartic = quartic5
#Kr_list0 = construct_Kr(quartic)
#Kr_list = []
#for Kr in Kr_list0:
#    if Kr not in Kr_list:
#        Kr_list.append(Kr)
#print(Kr_list)

Kr_list = Kr_d10

Kr_list_short = []
for i in range(0,100):
    Kr_list_short.append(Kr_list[i])
Kr_list = Kr_list_short

print(datetime.now())

K_list = CM_one_sextic_from_Kr(Kr_list)

print(datetime.now())

print(K_list)