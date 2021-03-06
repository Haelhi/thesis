lst_f = [i.strip('\n').split() for i in open('cubics.txt')]
lst_f.pop(); lst_f.pop()

cubics = []
for c in lst_f:
    entries = [int(x) for x in c[2][1:-1].split(',') if x != '']
    entries2 = [int(y) for y in c[1][1:-1].split(',') if x != '']
    cubics.append([entries[::-1],entries2])

def check_hessian(hess):
    for i in hess:
        if hess[0] != i and hess[0] != -i:
            return False
    return True
    
Kplus_C3 = [i for i in cubics if check_hessian(i[1])==True]
Kplus_S3 = [i for i in cubics if check_hessian(i[1])==False]

o = open('cubics_7-8_C3.sage','a')
o.write('listC3 = ')
o.write(str(Kplus_C3))
o.close()
o = open('cubics_7-8_S3.sage','a')
o.write('listS3 = ')
o.write(str(Kplus_S3))
o.close()