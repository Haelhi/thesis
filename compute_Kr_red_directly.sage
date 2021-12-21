load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
from datetime import datetime 

load('functions.sage')
load('Data/quartics10.sage')

pari.allocatemem(83613065216)

quartic = quartic10

construct_Kr_red_write(quartic)