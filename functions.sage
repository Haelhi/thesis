
# ----------------------- NEW METHOD FOR ALL QUADRATICS SUCH THAT dk|dKrplus -------------------------

# Function construct_Kr: computes composite field Kr = kKrplus with all k such that dk|dKrplus
# INPUT: LMFDB output list of quartic fields
# OUTPUT: list with elements [polynomial of Kr, polynomial of k, class number of Krplus]
def construct_Kr(quartic):
    Kr_list = []
    for quart in quartic:
        Krplus.<o> = NumberField(quart)
        R.<t> = PolynomialRing(Krplus)
        dKrplus = Krplus.discriminant()
        kprimes = list(dKrplus.factor())
        for i in range(1,len(kprimes)+1):
            if i == 1:
                prime_combos = kprimes
                for p in prime_combos:
                    k.<u> = QuadraticField(-p[0])
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
                        dKr = Kr.discriminant()
                        if dKr == dKrplus^2:
                            Kr_list.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
            else:
                prime_combos = list(combinations(kprimes,i))
                for p in prime_combos:
                    d = 1
                    for n in range(0,len(p)):
                        d = d * p[n][0]
                    k.<u> = QuadraticField(-d)
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
                        dKr = Kr.discriminant()
                        if dKr == dKrplus^2:
                            Kr_list.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
    return Kr_list

# Can this be made more efficient with isomorphism check?
# Function compute_K_from_Kr_sub computes sextic fields K from reflex fields Kr up to isomorphism
# INPUT: list of non-normal octic reflex fields Kr with non-normal subfield Krplus
# OUTPUT: list with [Kr, Phir, K, k, class number of Krplus]
def compute_K_from_Kr_sub(list_Kr):
    list_sextic = []
    check_isom = 0
    for Kr in list_Kr:
        Krcm = CM_Field(Kr[0])
        k = CM_Field(Kr[1])
        Phir_set = Krcm.CM_types()
        for Phir in Phir_set:
            K = Phir.reflex_field()
            if K.g() == 3:
                if len(list_sextic) == 0:
                    list_sextic.append([Krcm,Phir,K,k,Kr[2]])
                else:
                    N.<a> = NumberField(K.polynomial())
                    for m in list_sextic:
                        n.<b> = NumberField(m[0].polynomial())
                        if N.is_isomorphic(n) == True:
                            check_isom = 1
                            break
                    if check_isom == 0:
                        list_sextic.append([Krcm,Phir,K,k,Kr[2]])
                    else:
                        check_isom = 0
    return list_sextic


# ----------------------- CHECK IF SEXTIC FIELD IS OF CORRECT FORM -------------------------
def check_list_K(list_sextic):
    # Check if sextic fields with polynomial in list_sextic are of type I:
    # not Galois, no imaginary quadratic subfield and Galois maximally totally real subfield
    list_K = []
    for K in list_sextic:
        Knf.<a> = NumberField(K[0].polynomial())
        if Knf.is_galois() == False:
            if all(k[0].degree()!=2 for k in Knf.subfields()):
                for Kplus in Knf.subfields():
                    if Kplus[0].degree() == 3:
                        if Kplus[0].is_galois() == true:
                            if len(Kplus[0].real_embeddings()) == 3:
                                list_K.append(K) 
    return list_K

# ----------------------- FIRST CHECK CM ONE THEN COMPUTE K ---------------------------

def CM_one_sextic_from_Kr(Kr_list):
    K_list = []
    for Kr_pol in Kr_list:
        print(Kr_pol[0])
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
                if K.g() == 1:
                    break
            if i == 6:
                break
    return K_list

    

            
    
# ----------------------- CHECK IF FIELD ISOMORPHIC TO OTHER FIELDS IN LIST ------------------------
# Check if list contains isomorphic K
# INPUT: (CM-field K, list of polynomials of computed fields K)
# OUTPUT: "True" if there are no isomorphic K in the existing list, "False" if there are
def check_isomorphic_K_in_list(K,list_sextic):
    check_isom = 0
    if len(list_sextic) == 0:
        return "True"
    else:
        N.<a> = NumberField(K.polynomial())
        for k in list_sextic:
            n.<b> = NumberField(k[0].polynomial())
            if N.is_isomorphic(n) == True:
                check_isom = 1
                break
    if check_isom == 0:
        return "True"
    return "False"
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

# ----------------------- WORK AROUND REPRESENTATIVE_PRIME BUG ------------------------
def execute_rep_prime(nb):
    try:
        return representative_prime(norm_bound=nb)
    except:
        nb = nb*10
        return execute_rep_prime(nb)