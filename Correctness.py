from sage.schemes.hyperelliptic_curves import jacobian_morphism
from sage.all import random_prime
import random
from sage.rings.finite_rings.finite_field_constructor import FiniteFieldFactory 
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve


def RandomCurve_RAM(g,q):
    
    GF = FiniteFieldFactory(q)
    x = GF(q)['x'].gen()
   
    while True:
        f = x**(2*g+1)
        h = x**(g +1)
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
    
    print("Generating Full Curves")
    genus = g
    testingCurvesFull = []
    
    i = 0
    while i < num:
        H,f,h = RandomCurve_RAM(g,q)
        #print(f,h)
        
        J = H.jacobian()
        J = J(J.base_ring())
        
        div1 = 0
        div2 = 0
        
        divs = []
        
        j = 0
        k = 0
        while j < g:
            try:
                divs.append(J(H.lift_x(k)))
                j += 1
            except:
                pass
            k += 1
        
        div1 = divs[0]
        for d in divs[1:]:
            div1 += d
        
        divs = []
        j = 0
        while j < g:
            try:
                divs.append(J(H.lift_x(k)))
                j += 1
            except:
                pass
            k += 1
        
        div2 = divs[0]
        for d in divs[1:]:
            div2 += d
        
        i += 1
    
        
        testingCurvesFull.append((H,f,h,div1,div2))
        
    
    print("Testing add_NUCOMP")
    i = 0
    #print(div1,div2)
    for curve in testingCurvesFull:
        H = curve[0]
        div1 = curve[3]
        div2 = curve[4]
        X = div1.parent()
        
        f, h = H.hyperelliptic_polynomials()
        genus = H.genus()
        
        #if (testingCurvesFull[i][3] != testingCurvesFull[i][3].booger): print('false')
        u,v,w = jacobian_morphism.add_NUCOMP(div1,div2,f,h,genus,div1.w,div2.w)
        temp1 = (u,v)
        D1 = jacobian_morphism.cantor_composition(div1,div2,f,h,genus)
        if D1[0].degree() > genus:
            D1 = jacobian_morphism.cantor_reduction(D1[0], D1[1], f,h,genus)
        
        if (temp1 != D1): print("add_NUCOMP - Failed",temp1,D1)
        
        i += 1
    
    
  
    
    
    print("Testing NUDUPL")
    i = 0
    for curve in testingCurvesFull:
        H = curve[0]
        div1 = curve[3]
        div2 = curve[4]
        X = div1.parent()
        
        f, h = H.hyperelliptic_polynomials()
        genus = H.genus()
        
        #if (testingCurvesFull[i][3] != testingCurvesFull[i][3].booger): print('false')
        u,v,w = jacobian_morphism.add_NUCOMP(div1,div1,f,h,genus,div1.w,div1.w)
        temp1 = (u,v)
        D1 = jacobian_morphism.cantor_composition(div1,div1,f,h,genus)
        if D1[0].degree() > genus:
            D1 = jacobian_morphism.cantor_reduction(D1[0], D1[1], f,h,genus)
        
        if (temp1 != D1): print("add_NUCOMP - Failed",temp1,D1)
        
        i += 1
    
    
    
    
    
    print("Generating Simple Curves")
    
    testingCurvesSimple = []
    genus = g
    testingCurvesFull = []
    
    i = 0
    while i < num:
        H,f = RandomCurve_RAM_simple(g,q)
        #print(f,h)
        
        J = H.jacobian()
        J = J(J.base_ring())
        
        div1 = 0
        div2 = 0
        
        divs = []
        
        j = 0
        k = 0
        while j < g:
            try:
                divs.append(J(H.lift_x(k)))
                j += 1
            except:
                pass
            k += 1
        
        div1 = divs[0]
        for d in divs[1:]:
            div1 += d
        
        divs = []
        j = 0
        while j < g:
            try:
                divs.append(J(H.lift_x(k)))
                j += 1
            except:
                pass
            k += 1
        
        div2 = divs[0]
        for d in divs[1:]:
            div2 += d
        
        i += 1
    
        
        testingCurvesFull.append((H,f,div1,div2))
    
    #print(div1,div2)
    print("Testing add_NUCOMP_simple")
    for curve in testingCurvesSimple:
        H = curve[0]
        div1 = curve[2]
        div2 = curve[3]
        X = div1.parent()
        
        f, h = H.hyperelliptic_polynomials()
        if(h != 0): print(h)
        genus = H.genus()
        #if (testingCurvesFull[i][3] != testingCurvesFull[i][3].booger): print('false')
        u,v,w = jacobian_morphism.add_NUCOMP_simple(div1,div2,f,genus,div1.w,div2.w)
        temp1 = (u,v)
        
        #temp1 = jacobian_morphism.JacobianMorphism_divisor_class_field(X,(u,v),w = w, check=False)
        
        D1 = jacobian_morphism.cantor_composition_simple(div1,div2,f,genus)
        if D1[0].degree() > genus:
            D1 = jacobian_morphism.cantor_reduction_simple(D1[0], D1[1], f,genus)
            
        #D1 = jacobian_morphism.JacobianMorphism_divisor_class_field(X,D1, check=False)    
        
        if (temp1 != D1): print("add_NUCOMP_simple - Failed",temp1,D1)
        
    print("Testing add_NUDUPL_simple")
    for curve in testingCurvesSimple:
        H = curve[0]
        div1 = curve[2]
        
        X = div1.parent()
        f, h = H.hyperelliptic_polynomials()
        if(h != 0): print(h)
        genus = H.genus()
        u,v,w = jacobian_morphism.add_NUCOMP_simple(div1,div1,f,genus,div1.w,div1.w)
        temp1 = (u,v)
        
        D1 = jacobian_morphism.cantor_composition_simple(div1,div1,f,genus)
        if D1[0].degree() > genus:
            D1 = jacobian_morphism.cantor_reduction_simple(D1[0], D1[1], f,genus)
        
        if (temp1 != D1): print("add_NUDUPL_simple - Failed",temp1,D1)
    
    
#Black Box Testing    
fieldsizes = [16]
numberOfTests = 1000


for fieldsize in fieldsizes:
    print("Field Size: {0}".format(fieldsize))
    prime = random_prime(2**fieldsize)
    #prime = 127
    
    for g in range(1,31):
        print("Genus: {0}".format(g))
        #print(prime)
        runTest(numberOfTests,g,prime)
        

    