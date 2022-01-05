
# There are 3 functions in this file, they eliminate fields wrt the occurance of different types of decompositions


# We can use the following function BEFORE we construct Kr. The function is only meaningful when p divides BOTH dk and dKrplus

def elimination_wrt_decomposition_1_3_7_and_15(p, Krplus_min_poly, dKrplus, dk): 
    """
    
    Tests whether Kr/Q has a decomposition type 1, 3, 7 or 15 as in Table ??
    # This only happens if p divides dk and dKrplus
    
    This is Step ?? in Algorithm ??

    INPUT:
    ''Krplus_min_poly'' -- minimal polynomial of a totally real quartic
    ''dk''  -- discriminant of k
    ''dKrplus'' -- discriminant of Krplus
    ''p'' -- prime number dividing both dk and dKrplus

    OUTPUT:

    - Returns True if p does not have decomposition type 1, 3, 7 or 15 as in Table ??

    """
    
    if dk % p != 0 or dKrplus % p != 0:
        print(p, dk.factor(), dKrplus.factor())
        raise RuntimeError("This is not a correct prime since p divides dk and dKrplus")
    
    T.<x> = PolynomialRing(GF(p))
    poly_over_finite = T(Krplus_min_poly)
    # Eliminates fields with primes of type 1, 3, or 7:
    if len(poly_over_finite.factor()) == 1 and poly_over_finite.factor()[0][1] == 4:
        if len(Krplus.factor(p)) == 1 and Krplus.factor(p)[0][1] == 4:
            return False
    # Eliminates fields with primes of type 15:
    if len(poly_over_finite.factor()) == 2:
        if (poly_over_finite.factor()[0][1] == 1 and poly_over_finite.factor()[1][1] == 3) or (poly_over_finite.factor()[0][1] == 3 and poly_over_finite.factor()[1][1] == 1):
            if len(Krplus.factor(p)) == 2:
                if (Krplus.factor(p)[0][1] == 1 and Krplus.factor(p)[1][1] == 3) or (Krplus.factor(p)[0][1] == 3 and Krplus.factor(p)[1][1] == 1):
                    return False
    else:
        pass
    return True
    

# We can use the following function AFTER we construct Kr. Prime p divides dKrplus but NOT dk. (You dont 
# Unfortunately we don't eliminate much :/ 

def elimination_wrt_decomposition_5_6_and_18(p, Kr_min_poly, dKrplus, dk): 
    """
    
    Tests whether Kr/Q has a decomposition type 5, 6 or 18 as in Table ??
    # This only happens if p divides dKrplus and does not divide dk
    
    This is Step ?? in Algorithm ??

    INPUT:
    ''Kr_min_poly'' -- minimal polynomial of an octic CM field
    ''dk''  -- discriminant of k
    ''dKrplus'' -- discriminant of Krplus
    '' p '' -- prime number dividing dKrplus but not dk

    OUTPUT:

    - Returns True if there are no primes in Kr/Q with decomposition type 5, 6 or 18 as in Table ??

    """
    
    if not (dk % p != 0 and dKrplus % p == 0):
        print(p, dk.factor(), dKrplus.factor())
        raise RuntimeError("This is not a correct prime since p divides dKrplus and does not divide dk")
    
    T.<x> = PolynomialRing(GF(p))
    Kr_poly_over_finite = T(Kr_min_poly)
    # Eliminates fields with primes of type 5, 6 or 18:
    if len(Kr_poly_over_finite.factor()) == 2 and Kr_poly_over_finite.factor()[0][1] == 4:
        if len(Kr.factor(p)) == 2 and Kr.factor(p)[0][1] == 4:
            return False
    else:
        pass
    return True
	
	
# We can use the following function AFTER we construct Kr. Prime p divides BOTH dKrplus and dk.

def elimination_wrt_decomposition_19(p, Krplus_min_poly, Kr_min_poly, dKrplus, dk): 
    """
    
    Tests whether Kr/Q has a decomposition type 19 as in Table ??
    # This only happens if p divides dKrplus and does not divide dk
    
    This is Step ?? in Algorithm ??

    INPUT:
    ''Kr_min_poly'' -- minimal polynomial of an octic CM field
    ''Krplus_min_poly'' -- minimal polynomial of the totally real field in Kr
    ''dk''  -- discriminant of k
    ''dKrplus'' -- discriminant of Krplus
    '' p '' -- prime number dividing dKrplus but not dk

    OUTPUT:

    - Returns True if there are primes in Kr/Q with decomposition type 19 as in Table ??

    """
    
    if not (dk % p == 0 and dKrplus % p == 0):
        print(p, dk.factor(), dKrplus.factor())
        raise RuntimeError("This is not a correct prime since p divides dKrplus and does not divide dk")
    
    T.<x> = PolynomialRing(GF(p))
    Kr_poly_over_finite = T(Kr_min_poly)
    Krplus_poly_over_finite = T(Krplus_min_poly)
    # Eliminates fields with primes of type 19:
    if len(Krplus_poly_over_finite.factor()) == 2 and Kr_poly_over_finite.factor()[0][1] == 4:
        if len(Krplus.factor(p)) == 2 and Kr.factor(p)[0][1] == 4:
            return False
    else:
        pass
    return True