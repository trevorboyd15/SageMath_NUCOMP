from sage.schemes.hyperelliptic_curves import jacobian_morphism
from sage.all import random_prime
import random
import time
from sage.rings.finite_rings.finite_field_constructor import FiniteFieldFactory 
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.all import timeit




def RandomCurve_RAM(g,q):
    
    GF = FiniteFieldFactory(q)
    x = GF(q)['x'].gen()
   
    while True:
        f = x**(2*g+1)
        h = x**(g + 1)
        for i in range(0,g + 1):
            h += random.randrange(0,q)*x**i
        
        for i in range(0,2*g+1):
            f += random.randrange(0,q)*x**i
        try:
            H = HyperellipticCurve(f,h)
            if H.genus() == g:
                
                return H,f,h
        except:
            pass
                   
def RandomCurve_RAM_simple(g,q):
    
    GF = FiniteFieldFactory(q)
    x = GF(q)['x'].gen()
   
    while True:
        f = x**(2*g+1)
        
        for i in range(0,2*g+1):
            f += random.randrange(0,q)*x**i
        try:
            
            H = HyperellipticCurve(f)
            if H.genus() == g:
                return H,f
        except:
            pass

def runTest(num,g,q):

    testingCurvesFull = []
    H,f,h = RandomCurve_RAM(g,q)
        
    J = H.jacobian()
    J = J(J.base_ring())
        
    div1 = 0
    div2 = 0
    divs = []
        
    i = 0 
    j = 0
    
    #Generate Divisor 
    while i < g:
        try:
            divs.append(J(H.lift_x(j)))
            i += 1
        except:
            pass
        j += 1
    div1 = divs[0]
        
        
    for d in divs[1:]:
        div1 += d
            
            
    divs = []
    
    #Generate Second Divisor
    i = 0 
    while i < g:
        try:
            divs.append(J(H.lift_x(j)))
            i += 1
        except:
            pass
        j += 1
    div2 = divs[0]
    for d in divs[1:]:
        div2 += d
            
            
        
    testingCurvesFull.append((H,f,h,div1,div2))
    
    
    print("Testing Cantor_Full")
    
    X = div1.parent()
    H = testingCurvesFull[0][0]
    f, h = H.hyperelliptic_polynomials()
    genus = H.genus()
    
    t = time.process_time()
    while i < num:
        for curve in testingCurvesFull:
           
            D1 = jacobian_morphism.cantor_composition(div1,div2,f,h,genus)
            if D1[0].degree() > genus:
                D1 = jacobian_morphism.cantor_reduction(D1[0], D1[1], f,h,genus)
            
            
        i += 1
    t2 = time.process_time()
    Cantor_Full_time = t2 - t
    
    print("Testing add_NUCOMP")
    i = 0
    t = time.process_time()
    while i < num:
        u,v,w = jacobian_morphism.add_NUCOMP(div1,div2,f,h,genus,div1.w,div2.w)
        i += 1
    t2 = time.process_time()
    NUCOMP_Full_time = t2 - t
    
    
    print("Testing NUDUPL")
    X = div1.parent()
    H = testingCurvesFull[0][0]
    f, h = H.hyperelliptic_polynomials()
    genus = H.genus()
    
    
    i = 0
    t = time.process_time()
    while i < num:
        u,v,w = jacobian_morphism.add_NUCOMP(div1,div1,f,h,genus,div1.w,div1.w)
        i += 1
    t2 = time.process_time()
    
    NUDUPL_time = t2-t
    
    
    print("Generating Simple Curves")
    testingCurvesSimple = []
    
    i = 0
    
    H,f = RandomCurve_RAM_simple(g,q)
        
    J = H.jacobian()
    J = J(J.base_ring())
        
    div1 = 0
    div2 = 0
    divs = []
        
    i = 0 
    j = 0
    while i < g:
        try:
            divs.append(J(H.lift_x(j)))
            i += 1
        except:
            pass
        j += 1
    div1 = divs[0]
           
    for d in divs[1:]:
        div1 += d
  
    divs = []
        
    i = 0 
    while i < g:
        try:
            divs.append(J(H.lift_x(j)))
            i += 1
        except:
            pass
        j += 1
    div2 = divs[0]
    for d in divs[1:]:
        div2 += d
            
    testingCurvesSimple.append((H,f,div1,div2))
        
        
    print("Testing NUCOMP_Simple")
    i = 0
    t = time.process_time()
    while i < num:
        u,v,w = jacobian_morphism.add_NUCOMP_simple(div1,div2,f,genus,div1.w,div2.w)
        #D1 = jacobian_morphism.JacobianMorphism_divisor_class_field(X,(u,v),w = w, check=False)
        i += 1
    t2 = time.process_time()
    
    NUCOMP_Simple_time = t2-t
        
    
    print("Testing Cantor_Simple")    
    i = 0
    t = time.process_time()
    while i < num:
        D1 = jacobian_morphism.cantor_composition_simple(div1,div2,f,genus)
        if D1[0].degree() > genus:
            D1 = jacobian_morphism.cantor_reduction_simple(D1[0], D1[1], f,genus)
        i += 1
    t2 = time.process_time()
    Cantor_Simple_time = t2 - t

    print("Testing NUDUPL_simple")
    i = 0
    t = time.process_time()
    while i < num:
        u,v,w = jacobian_morphism.add_NUCOMP_simple(div1,div1,f,genus,div1.w,div1.w)
        #D1 = jacobian_morphism.JacobianMorphism_divisor_class_field(X,(u,v),w = w, check=False)
        i += 1
    t2 = time.process_time()
    NUDUPL_simple_time = t2-t
    
    return NUCOMP_Full_time,NUDUPL_time,NUCOMP_Simple_time,NUDUPL_simple_time,Cantor_Full_time,Cantor_Simple_time
    

fieldsizes = [8,16,32,64,128]
numberOfTests = 1000


for fieldsize in fieldsizes:
    print("Field Size: {0}".format(fieldsize))
    prime = random_prime(2**fieldsize)
    
    f = open("res{0}.csv".format(fieldsize),'w')
    f.write("Genus,NUCOMP_Full_time,NUDUPL_time,NUCOMP_Simple_time,NUDUPL_simple_time,Cantor_Full_time,Cantor_Simple_time\n")
    for g in range(1,51):
        print("Genus: {0}".format(g))
        NUCOMP_Full_time,NUDUPL_time,NUCOMP_Simple_time,NUDUPL_simple_time,Cantor_Full_time,Cantor_Simple_time = runTest(numberOfTests,g,prime)
        f.write("{0},{1},{2},{3},{4},{5},{6}\n".format(g,NUCOMP_Full_time,NUDUPL_time,NUCOMP_Simple_time,NUDUPL_simple_time,Cantor_Full_time,Cantor_Simple_time))
        f.flush()
    f.close()
        
    
    
    
    
    
    