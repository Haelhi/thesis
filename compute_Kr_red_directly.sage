load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")
from datetime import datetime 

load('functions.sage')
load('Data/quartics10.sage')

pari.allocatemem(83613065216)

quartic = quartic10

(Kr_one,Kr_prime) = construct_Kr_reduced(quartic)

o = open('Data/Kr_d10_p_div_dKplus_red.sage', 'a')
o.write('Kr_d10_one = [')
o.write(str(Kr_one))
o.write(']')
o.write('\n\n\n')
o.write('Kr_d10_prime = [')
o.write(str(Kr_prime))
o.write(']')
o.close()