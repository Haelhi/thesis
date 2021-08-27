load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations
from datetime import datetime

pari.allocatemem(80000000000)

load('functions.sage')
load('data_K.sage')

list_K = add_K
comp_list = []

for K in list_K:
    Kcm = CM_Field(K)
    Phi = Kcm.CM_types()[0]
    Krcm = Phi.reflex_field()
    comp_list.append([K,Krcm.polynomial()])

print(comp_list)