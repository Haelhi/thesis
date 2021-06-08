load("https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage")

from itertools import combinations

load('functions.sage')
load('data_subfields.sage')


# quartic = quartic_2
# Kr_list0 = construct_Kr(quartic)
# print(len(Kr_list0))
# Kr_list = []
# for Kr in Kr_list0:
#     if Kr not in Kr_list:
#         Kr_list.append(Kr)
# print(len(Kr_list))

# print(f'[{Kr_list[0]}', *Kr_list[1:-1], f'{Kr_list[-1]}]', sep = ',\n')

# print(len(Kr_list))

# Kr_list = list_Kr_1

# Kr_list_short = []
# for i in range(0,10):
#     Kr_list_short.append(list_Kr_1[i])
# Kr_list = Kr_list_short

Kr_list = list_Kr_2
K_list = compute_K_from_Kr_sub(Kr_list)
K_CM_one = K_cm_clno_one_list(K_list)
print(K_CM_one)
# print(f'[{K_CM_one[0]}', *K_CM_one[1:-1], f'{K_CM_one[-1]}]', sep = ',\n')