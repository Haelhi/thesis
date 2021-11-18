
# ----------------------- CONSTRUCT Kr FROM QUARTIC FIELDS -------------------------

# CASE 1: dKr/dKrplus^2 == 1
# Function construct_Kr: computes composite field Kr = kKrplus with all k such that dk|dKrplus 
# INPUT: LMFDB output list of quartic fields
# OUTPUT: list with elements [polynomial of Kr, polynomial of k, class number of Krplus]
def construct_Kr_1(quartic):
    Kr_list = []
    for poly in quartic:
        Krplus.<a> = NumberField(poly)
        R.<t> = PolynomialRing(Krplus)
        dKrplus = Krplus.discriminant()
        primes = list(dKrplus.factor())
        for i in range(1,len(primes)+1):
            if i == 1:
                p_combos = [[p] for p in primes]
            else:
                p_combos = list(combinations(primes,i))
            for p in p_combos:
                d = 1
                for n in range(0,len(p)):
                    d = d * p[n][0]
                k.<b> = QuadraticField(-d)
                dk = k.discriminant()
                if dKrplus % dk == 0:
                    kKrplus.<y> = Krplus.extension(k.polynomial())
                    Kr.<z> = kKrplus.absolute_field()
                    dKr = Kr.discriminant()
                    if dKr/dKrplus^2 == 1:
                        Kr_list.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
    return Kr_list


# CASE 2: dKr/dKrplus^2 == p and just use k = QQ(sqrt(-p))
# Function construct_Kr: computes composite field Kr = kKrplus with all k such that dk|dKrplus 
# INPUT: LMFDB output list of quartic fields
# OUTPUT: list with elements [polynomial of Kr, polynomial of k, class number of Krplus]

def construct_Kr_complete(quartic):
    Kr_list_one = []
    Kr_list_prime = []
    for poly in quartic:
        Krplus.<a> = NumberField(poly)
        R.<t> = PolynomialRing(Krplus)
        dKrplus = Krplus.discriminant()
        primes = list(dKrplus.factor())
        for P in primes:
            if P[0] % 4 == 3 or P[0] == 2:
                k.<b> = QuadraticField(-P[0])
                dk = k.discriminant()
                if dKrplus % dk == 0:
                    kKrplus.<y> = Krplus.extension(k.polynomial())
                    Kr.<z> = kKrplus.absolute_field()
                    dKr = Kr.discriminant()
                    if dKr/dKrplus^2 == 1:
                            Kr_list_one.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
                    if dKr/dKrplus^2 != 1:
                            Kr_list_prime.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
    return (Kr_list_one, Kr_list_prime)

# ----------------------- PARALLEL: COMPUTE K AND WRITE IN FILE ---------------------------

def CM_one_sextic_from_Kr_write(Kr_list):
    j = 0
    K_list = []
    for Kr_pol in Kr_list:
        Kr = CM_Field(Kr_pol[0])
        clKr = Kr.class_group()
        gens = clKr.gens()
        rep_gens = []
        for I in gens:
            M = Kr.minkowski_bound().numerical_approx().ceil()
            Ip = I.representative_prime(norm_bound=M)
            rep_gens.append(Ip)
        Phir_set = Kr.CM_types(equivalence_classes=True)
        i = 0
        for Phir in Phir_set:
            i = i + 1
            if test_CM_cl_nr_one_with_class_group(rep_gens, Phir) == True:
                K = Phir.reflex_field()
                if K.g() == 3:
                    polyK = pari.polredabs(K.polynomial())
                    o = open('parallel_OUTPUT.sage','a') 
                    o.write(str([polyK,Kr_pol[0],Kr_pol[1],Kr_pol[2]]))
                    o.write(', ')
                    o.close()
                    break
            if i == 6:
                break
        j = j + 1
    return K_list


# ----------------------- FIRST CHECK CM ONE THEN COMPUTE K ---------------------------

def CM_one_sextic_from_Kr(Kr_list):
    j = 0
    K_list = []
    for Kr_pol in Kr_list:
        print(Kr_pol[0])
        print(j)
        Kr = CM_Field(Kr_pol[0])
        clKr = Kr.class_group()
        gens = clKr.gens()
        rep_gens = []
        for I in gens:
            M = Kr.minkowski_bound().numerical_approx().ceil()
            Ip = I.representative_prime(norm_bound=M)
            rep_gens.append(Ip)
        Phir_set = Kr.CM_types(equivalence_classes=True)
        i = 0
        for Phir in Phir_set:
            print("Phir")
            i = i + 1
            if test_CM_cl_nr_one_with_class_group(rep_gens, Phir) == True:
                K = Phir.reflex_field()
                print(K.polynomial())
                if K.g() == 3:
                    K_list.append([pari.polredabs(K.polynomial()),Kr_pol[0],Kr_pol[1],Kr_pol[2]])
                    break
            if i == 6:
                break
        j = j + 1
    return K_list

# ----------------------- PARALLEL COMPUTATION FUNCTIONS ---------------------------

# INPUT: L: list to be split up; n: number of chunks
# OUTPUT: list of chunks from L
def divide_into_chunks(L,n):
    Chunk_list = []
    size = len(L)
    chunk_size = ceil(size/n)
    s = [L[j * chunk_size : (j+1)*chunk_size] for j in range(n)]
    for r in s:
        if r != []:
            Chunk_list.append(r)
    return Chunk_list

@parallel(4)
def CM_one_from_partial_Kr(Kr_list,x):
    Kr_list_short = []
    for i in range(x[0],x[1]):
        Kr_list_short.append(Kr_list[i])
    K_list = CM_one_sextic_from_Kr(Kr_list_short)
    return K_list
    
    
# ----------------------- CM CLASS NUMBER ONE FUNCTION WITH BACH BOUND -------------------------

def test_CM_cl_nr_one_with_bach_bound(Kr, Phir):
    """
    Test under GRH whether K has CM class number one.
    The output True means that either K has CM class number one or GRH is false.
    The output False means that K does not have CM class number one.

    Uses splitting of small primes as well as Weil Q-number enumeration.

    This is Step 3 in Algorithm 3. 
    
    INPUT:

    - ''Kr'' - a octic non-Galois CM field  
    - ''Phir'' - a CM type of Kr

    OUTPUT:
    
     - See above.

    """
    bach_bound = floor(Kr.bach_bound()) 
    for l in prime_range(bach_bound + 1): 
        for l_1 in Kr.factor(l):
            norm_l_1 = l_1[0].norm()
            if norm_l_1 <= bach_bound:
                if a_to_mu(Phir, l_1[0]) is None: 
                    return False 
    return True
    
    
# ----------------------- CM CLASS NUMBER ONE FUNCTION -------------------------

def test_CM_cl_nr_one_with_class_group_2(Kr, Phir):
    clKr = Kr.class_group()
    gens = clKr.gens()
    for I in gens:
        M = Kr.minkowski_bound().numerical_approx().ceil()
        Ip = I.representative_prime(norm_bound=M)
        if a_to_mu(Phir, Ip) is None: 
            return False 
    return True

# Function test_CM_cl_nr_one_with_class_group: test if K is a CM-class number one field
# INPUT: (prime representatives of generators of Cl_Kr, CM-type Phir of Kr)
# OUTPUT: True if CM-class number one field.
def test_CM_cl_nr_one_with_class_group(rep_gens, Phir):
    """
    Checks whether the given non-biquadratic quartic CM field is a PQ1-field.
    (In other words, checks whether I_0(Phir)/P_{Kr} = Cl_{Kr}.)

    This function uses the class group of K, hence is only fast for fields K of small discriminant.
    We use it only after first convincing ourselves in other ways whether that the answer will be 'yes'.

    This is Step 3 in Algorithm 3. # TODO: check later whether these references are correct.

    INPUT: 

    '' Kr '' - non-biquadratic quartic CM field
    '' Phir '' - is a CM type of Kr
    '' K '' - the reflex field of (Kr ,Phir)

    OUTPUT: returns True if K is a PQ1 field
    
    METHOD: 
    
    Step 1. Compute the class group of $Kr$
    Step 2. For each prime ideal representatives "Ip" of the classes compute N_{Phir}(Ip) and check if it is
    generated by a Weil N(Ip)-number as follows.
    a. First check whether it is a principal ideal.
    b. If yes, say N_{Phir}(Ip) = (alpha0), then v = compute alpha0*alpha0bar / N(Ip).
    c. Now we need to check whether the totally positive unit v is a relative norm from K to K_plus.
       We can easily do this using a_to_mu from the ReCip package,
       but as O_K^* = mu_K*O_{K_plus}^*, [13, Lemma 1] # TODO: update reference later as [13] can get a different number when we make changes.
       it is also equivalent to being a square in K_plus, which is equivalent to being a square in K.
       So we test whether v is a square.
    """
    
#     clKr = Kr.class_group()
#     gens = clKr.gens()
#     for I in gens:
#         M = Kr.minkowski_bound().numerical_approx().ceil()
#         Ip = I.representative_prime(norm_bound=M)
    for Ip in rep_gens:
        if a_to_mu(Phir, Ip) is None: 
            return False 
    return True
    
def test_CM_cl_nr_one_with_class_group_og(Kr, Phir):
    clKr = Kr.class_group()
    gens = clKr.gens()
    for I in gens:
        M = Kr.minkowski_bound().numerical_approx().ceil()
        Ip = I.representative_prime(norm_bound=M)
        if a_to_mu(Phir, Ip) is None: 
            return False 
    return True


# ----------------------- COMPUTE CM CLASS NUMBER ONE FIELDS K -------------------------
# INPUT: List of sextic reflex fields K of Kr
# OUTPUT: List of CM-class number one K [polynomial of K, hK, hk, hKr]
def K_cm_clno_one_list(K_list):
    list_CM_one = []
    for tup in K_list:
        Kr = tup[0]; Phir = tup[1]; K = tup[2]; k = tup[3]; hKrplus = tup[4]
        if test_CM_cl_nr_one_with_class_group(Kr,Phir) == True:
            hK = K.class_number()
            hKr = Kr.class_number()
            hKrstar = hKr/hKrplus
            hk = k.class_number()
            list_CM_one.append([K.polynomial(),hK,hk,hKrstar])
    return(list_CM_one)