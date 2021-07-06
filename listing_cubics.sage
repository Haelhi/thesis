f = [i.strip('\n').split() for i in open('cubics_8-9.txt')]
for elt in f:
    print(elt[2])