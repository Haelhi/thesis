# ----------------------- OLD METHOD USING CLASS NUMBER -------------------------
# Function Kr_from_quartic: computes composite fields of totally real non-normal quartic fields and
# imaginary quadratic fields of class number 1,2,4.
# INPUT: (list quartic fields, list imaginary quadratic fields)
# OUTPUT: list with elements [polynomial of Kr, dKr/dKrplus^2, hKrstar]
def Kr_from_quartic(quartic,iq):
    list_Kr = []
    check_isom = 0
    for quart in quartic:
        Krplus.<f0> = NumberField(quart[1])
        R0.<t0> = PolynomialRing(Krplus)
        hKrplus = Krplus.class_number()
        dKrplus = Krplus.discriminant()
        for quadr in iq:
            kKrplus.<f1> = Krplus.extension(quadr[1])
            Kr.<f> = kKrplus.absolute_field()
            dKr = Kr.discriminant()
            d = dKr/dKrplus^2
            if  d == 1:
                if len(list_Kr) == 0:
                    hKr = Kr.class_number()
                    hKrstar = hKr/hKrplus
                    list_Kr.append([Kr.absolute_polynomial(),d, hKrstar])
                else:
                    for m in list_Kr:
                        n.<b> = NumberField(m[0])
                        if Kr.is_isomorphic(n) == True:
                            check_isom = 1
                            break
                    if check_isom == 0:
                        hKr = Kr.class_number()
                        hKrstar = hKr/hKrplus
                        list_Kr.append([Kr.absolute_polynomial(),d, hKrstar])
                    else:
                        check_isom = 0
    return list_Kr

# Function compute_K_from_Kr computes sextic fields K from reflex fields Kr up to isomorphism
# INPUT: list of non-normal octic reflex fields Kr with non-normal subfield Krplus
# OUTPUT: list with [K,Kr,Phir] where K is the sextic reflex field of Kr and Phir the CM-type
def compute_K_from_Kr(list_Kr):
    list_sextic = []
    check_isom = 0
    for Kr in list_Kr:
        Krcm = CM_Field(Kr[0])
        Phir_set = Krcm.CM_types()
        for Phir in Phir_set:
            K = Phir.reflex_field()
            if K.g() == 3:
                if len(list_sextic) == 0:
                    list_sextic.append([K,Krcm,Phir])
                else:
                    N.<a> = NumberField(K.polynomial())
                    for k in list_sextic:
                        n.<b> = NumberField(k[0].polynomial())
                        if N.is_isomorphic(n) == True:
                            check_isom = 1
                            break
                    if check_isom == 0:
                        list_sextic.append([K,Krcm,Phir])
                    else:
                        check_isom = 0
    return list_sextic

# ----------------------- NEW METHOD FOR ALL QUADRATICS SUCH THAT dk|dKrplus -------------------------

# Function construct_Kr: computes composite field Kr = kKrplus with all k such that dk|dKrplus
# INPUT: LMFDB output list of quartic fields
# OUTPUT: list with elements [polynomial of Kr, polynomial of k, class number of Krplus]
def construct_Kr(quartic):
    Kr_list = []
    for quart in quartic:
        Krplus.<f0> = NumberField(quart[1])
        R0.<t0> = PolynomialRing(Krplus)
        dKrplus = Krplus.discriminant()
        kprimes = list(dKrplus.factor())
        for i in range(1,len(kprimes)+1):
            check_isom = 0
            if i == 1:
                prime_combos = kprimes
                for p in prime_combos:
                    k.<x> = QuadraticField(-p[0])
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
                        Kr_list.append([pari(Kr.polynomial()).polredabs(),pari(k.polynomial()).polredabs(),Krplus.class_number()])
            else:
                prime_combos = list(combinations(kprimes,i))
                for p in prime_combos:
                    d = 1
                    for n in range(1,i+1):
                        d = d * p[i-1][0]
                    k.<x> = QuadraticField(-d)
                    dk = k.discriminant()
                    if dKrplus % dk == 0:
                        kKrplus.<y> = Krplus.extension(k.polynomial())
                        Kr.<z> = kKrplus.absolute_field()
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
        Kr = CM_Field(Kr_pol[0])
        Phir_set = Kr.CM_types()
        for Phir in Phir_set:
            if test_CM_cl_nr_one_with_class_group(Kr, Phir) == True:
                K = Phir.reflex_field()
                if K.g() == 3:
                    K_list.append([K.polynomial(),Kr_pol[0],Kr_pol[1],Kr_pol[2]])
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
# ----------------------- CM CLASS NUMBER ONE FUNCTION -------------------------
# Function test_CM_cl_nr_one_with_class_group: test if K is a CM-class number one field
# INPUT: (CM-field K, CM-field Kr, CM-type Phir of Kr)
# OUTPUT: True if CM-class number one field.
def test_CM_cl_nr_one_with_class_group(Kr, Phir):
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
    
    clKr = Kr.class_group()
    gens = clKr.gens()
    for I in gens:
        Ip = I.representative_prime(norm_bound=50000000)
        # I = Phir.type_norm(Ip)
        # Now we need to check whether the ideal "I" is generated by a generator alpha such that alpha*alphabar is in QQ.
        # This can be done using
        #    if a_to_mu(Phir, I) is None: return False
        # from the ReCip package, which works for general CM fields.
        # But for primitive quartic CM fields, we give a direct check below.
        #if not I.is_principal():
        #    return False
        #Igen = I.gens_reduced()[0]
        #if not (Igen*Igen.conjugate() / Ip.absolute_norm()).is_square():
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