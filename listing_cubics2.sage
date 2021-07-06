lst_f = [i.strip('\n') for i in open('test_cubic.txt')]
lst_f.pop(); lst_f.pop()

list_entries = []
for c in lst_f:
    entries = [int(x) for x in c[1:-1].split(',') if x != '']
    list_entries.append(entries)
print(list_entries)