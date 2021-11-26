load('functions.sage')

def make_prime_list(start,end):
    q = start
    primes = []
    while q < end:
        q = next_prime(q)
        if q % 4 != 1:
            primes.append(q)
    return primes

def compute_max_hk(p_list):
    hk_list = []
    for p in p_list:
        k.<a> = QuadraticField(-p)
        hk = k.class_number()
        hk_list.append(hk)
    return(max(hk_list))

@parallel(30)
def list_of_max_hk(primes_chunks):
    max_hk = []
    for chunk in primes_chunks:
        hk_chunk = compute_max_hk([chunk])
        max_hk.append(hk_chunk)
    return(max_hk)

primes = make_prime_list(100000000,200000000)
primes_chunks = divide_into_chunks(primes,10000)
max_hk_list = list_of_max_hk(primes_chunks)
h = 1
for i in max_hk_list:
    if max(i[1]) > h:
        h = max(i[1])
print(h)