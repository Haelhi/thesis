f = [i.strip('\n').split() for i in open('cubics_6-7.txt')]
for elt in f:
    print(elt[2])