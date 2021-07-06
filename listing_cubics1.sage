lst_f = [i.strip('\n').split() for i in open('cubics.txt')]
lst_f.pop(); lst_f.pop()

list_entries = []
for c in lst_f:
    entries = [int(x) for x in c[2][1:-1].split(',') if x != '']
    list_entries.append(entries)
print(list_entries)