f = [i.strip('\n').split() for i in open('cubics_7-8.txt')]
for elt in f:
    print(elt[2])