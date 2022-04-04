from math import sqrt
import numpy as np

def unimodalni(tocka,h,f):

    left = tocka-h
    middle = tocka
    right = tocka+h
    step = 2

    fl = f(left)
    fm = f(middle)
    fr = f(right)

    fCalled = 3

    if fm<fr and fm<fl:
        return [left,right,fCalled]
    elif fm>fr:
        while True:

            left = middle
            middle = right
            fm = fr
            right = tocka + h*step
            step *= 2
            fr = f(right)
            fCalled += 1

            if not(fm > fr):
                return [left,right,fCalled]
    else:
        while True:
            right = middle
            middle = left
            fm = fl
            left = tocka - h*step
            step *= 2
            fl = f(left)
            fCalled += 1

            if not(fm>fl):
                return [left,right,fCalled]


def zlatniRez(f,e=10E-6,a=None,b=None,tocka=None,h=1,trace = False):

    if a is None and b is None and tocka is not None:
        a,b,fCalled = unimodalni(tocka,h,f)
    elif (a is None or b is None) and tocka is None:
        print("Pogreška prilikom poziva zlatnog reza")
        return

    k = 0.5*(sqrt(5)-1)

    korak = 1

    c = b - k * (b - a)
    d = a + k * (b - a)

    if trace == True:
        print("Korak 0: a = {}, c = {}, d = {}, b = {}".format(a,c,d,b))

    fc = f(c)
    fd = f(d)

    fCalled += 2

    while(b - a)>e:
        if fc<fd:
            b = d
            d = c
            c = b - k * (b - a)
            fd = fc
            fc = f(c)
            fCalled += 1
            
            if trace == True:
                print("Korak {}: a = {}, c = {}, d = {}, b = {}".format(korak,a,c,d,b))
                korak += 1

        else:
            a = c
            c = d
            d = a + k * (b - a) 
            fc = fd
            fd = f(d)
            fCalled += 1

            if trace == True:
                print("Korak {}: a = {}, c = {}, d = {}, b = {}".format(korak,a,c,d,b))
                korak += 1

    return [a,b,fCalled]

def euklidNorm(vektor):
    suma = 0
    for x in vektor:
        suma += x*x

    return sqrt(suma)
def gradientDescent(function,gradient,x0,precise = False, e = 10e-6):

    functionCalls = 0
    gradientCalculations = 0
    dimension = len(gradient)
    x = x0

    calculatedGradient = np.zeros(dimension)

    for index in range(0,dimension):
        calculatedGradient[index] = gradient[index](x)
    gradientCalculations += 1

    xValue = function(x)
    functionCalls += 1
    new_xValue = None
    diverging = 0


    while euklidNorm(calculatedGradient) > e:


        if not precise:
            x = x - calculatedGradient 
        else:
            left,right,fCall = zlatniRez(lambda lamb: function(x+lamb*calculatedGradient),e=e,tocka=0)
            lamb = (left+right)/2
            x = x + lamb*calculatedGradient
            functionCalls += fCall

        if euklidNorm(calculatedGradient) <= e:
            return functionCalls,gradientCalculations,x

        for index in range(0,dimension):
            calculatedGradient[index] = gradient[index](x)
        gradientCalculations += 1

        new_xValue = function(x)
        functionCalls += 1

        if new_xValue >= xValue:
            diverging += 1
        else:
            xValue = new_xValue
            diverging = 0

        if diverging >= 150:
            print("Postupak divergira")
            return functionCalls,gradientCalculations,x
        

    return functionCalls,gradientCalculations,x


def Newton_Raphson(function,gradient,H,x0,precise = False, e = 10e-6):
    
    functionCalls = 0
    gradientCalculations = 0
    hCalculations = 0
    dimension = len(gradient)
    rows = len(H)
    columns = len(H[0])
    x = x0

    calculatedGradient = np.zeros((dimension,1))
    calculatedH = np.zeros((rows,columns))

    for index in range(0,dimension):
        calculatedGradient[index][0] = gradient[index](x)
    gradientCalculations += 1

    for index1 in range(0,rows):
        for index2 in range(0,columns):
            calculatedH[index1][index2] = H[index1][index2](x)
    hCalculations += 1

    xValue = function(x)
    functionCalls += 1
    new_xValue = None
    diverging = 0


    while euklidNorm(calculatedGradient) > e:

        deltaX = np.linalg.inv(calculatedH).dot(calculatedGradient)

        if not precise:
            x = x - deltaX.T 
            x = x[0]

        else:
            deltaX = deltaX.T
            deltaX = deltaX[0]
            left,right,fCall = zlatniRez(lambda lamb: function(x+lamb*deltaX),e=e,tocka=0)
            lamb = (left+right)/2
            x = x + lamb*deltaX
            functionCalls += fCall

        if euklidNorm(calculatedGradient) <= e:
            return functionCalls,gradientCalculations,hCalculations,x

        for index in range(0,dimension):
            calculatedGradient[index][0] = gradient[index](x)
        gradientCalculations += 1

        for index1 in range(0,rows):
            for index2 in range(0,columns):
                calculatedH[index1][index2] = H[index1][index2](x)
        hCalculations += 1        

        new_xValue = function(x)
        functionCalls += 1

        if new_xValue >= xValue:
            diverging += 1
        else:
            xValue = new_xValue
            diverging = 0

        if diverging >= 100:
            print("Postupak divergira")
            return functionCalls,gradientCalculations,hCalculations,x
        

    return functionCalls,gradientCalculations,hCalculations,x


def boxovPostupak(function,eksplictina,implictina,dimension,x0,e=10e-6,alfa=1.3):

    for fun in eksplictina:
        if fun(x0) < 0:
            print("Početna točka ne zadovoljava ograničenja")
            return None
    for fun in implictina:
        if fun(x0) < 0:
            print("Početna točka ne zadovoljava ograničenja")
            return None

    skupTocaka = []
    skupTocaka.append(x0)

    xc = x0

    for index in range(0,2*dimension):
        tocka = np.random.rand(dimension)
        tocka = tocka*200
        tocka = tocka - 100
        iteracije = 0

        while True:
            


            toContinue = False
            for fun in implictina:
                if fun(tocka) < 0:
                    tocka = 0.5*(tocka+xc)
                    toContinue = True
                    break

            if toContinue:
                toContinue = False
                continue


            break

        skupTocaka.append(tocka)

        xc = sum(skupTocaka)/len(skupTocaka)


    while True:

        iteracije += 1

        if function(skupTocaka[0]) > function(skupTocaka[1]):
            h1 = 0
            h2 = 1
        else:
            h1 = 1
            h2 = 0

        for index in range(2,len(skupTocaka)):

            if function(skupTocaka[index]) > function(skupTocaka[h1]):
                h2 = h1
                h1 = index
            elif function(skupTocaka[index]) > function(skupTocaka[h2]):
                h2 = index

        xc = (sum(skupTocaka) - skupTocaka[h1])/(len(skupTocaka)-1)

        xr = (1+alfa)*xc-alfa*skupTocaka[h1]

        for index in range(0,dimension):

            if xr[index] < -100:
                xr[index] = -100
            elif xr[index] > 100:
                xr[index] = 100

        while True:

            toContinue = False
            for fun in implictina:
                if fun(xr) < 0:
                    xr = 0.5*(xr+xc)
                    toContinue = True
                    break

            if toContinue:
                toContinue = False
                continue

            break

        if function(xr) > function(skupTocaka[h2]):
            xr = 0.5*(xr+xc)

        skupTocaka.pop(h1)
        skupTocaka.append(xr)


        uvjetZaustavljanja = 0
        for index in range(0,len(skupTocaka)):
            uvjetZaustavljanja += pow((function(skupTocaka[index])-function(xc)),2)

        uvjetZaustavljanja /= len(skupTocaka)

        if uvjetZaustavljanja<e or iteracije >5000:

            maxf = function(skupTocaka[0])
            maxi = 0
            for index in range(1,len(skupTocaka)):
                value = function(skupTocaka[index])
                if value<maxf:
                    maxf = value
                    maxi = index

            return skupTocaka[maxi]
            

def istrazi(f, Xp, Dx):

    x = Xp

    #Iteriraj po komponentama vektora x
    for i in range(0,len(Xp)):
        P = f(x)

        #Povećaj i-tu komponentu za Dx
        x[i] = x[i] + Dx

        N = f(x)

        #Ako je vrijednost u novodobivenoj točki gora nego u početnoj, idi u drugom smjeru
        if N>P:

            x[i] = x[i] - 2*Dx

            N = f(x)
            
            #Ako je vrijednost i u drugom smjeru gora nego u početnoj, vrati točku na početnu
            if N>P:

                x[i] = x[i] + Dx

    return x

def HookJeevesUzOgranicenja(function,nejednakosti,jednakosti,x0,t=1,vx = 1,e = 10e-6,trace = False):
    Xp = np.array(x0).astype(float)
    Xb = np.array(x0).astype(float)
    Korak = 0
    

    while True:
        Korak += 1

        if Korak == 1:
            pass
        else:
            t = t+10


        f = lambda x: transformationFunction(function,nejednakosti,jednakosti,t,x)
        #Istrazi najbolju tocku
        Xn = istrazi(f,Xp.copy(),vx)
        Xn = np.array(Xn).astype(float)

        # Ako je novodobivena tocka bolja od prethodne, prihvati baznu tocku
        if f(Xn) < f(Xb):
            Xp = 2*Xn - Xb
            Xb = Xn
        #Ako novodobivena tocka nije bolja, smanji pomake(Dx)
        else:
            vx = vx/2
            Xp = Xb.copy()

        #Uvjet zaustavljanja
        if vx <= e:
            return Xb

def transformationFunction(function,nejednakosti,jednakosti,t,x):
    suma = 0

    suma += function(x)

    for nejed in nejednakosti:
        if nejed(x) < 0:
            return 9999999
        if nejed(x) != 0:
            suma -= (1/t)*np.log(nejed(x))

    for jed in jednakosti:
        suma += t*jed(x)**2

    return suma

def unutarnjaTocka(tocka,nejednakosti):
    suma = 0

    for nejed in nejednakosti:

        if nejed(tocka) < 0:
            suma -= nejed(tocka)

    return suma


def istrazi2(f, Xp, Dx):

    x = Xp

    #Iteriraj po komponentama vektora x
    for i in range(0,len(Xp)):
        P = f(x)

        #Povećaj i-tu komponentu za Dx
        x[i] = x[i] + Dx

        N = f(x)

        #Ako je vrijednost u novodobivenoj točki gora nego u početnoj, idi u drugom smjeru
        if N>P:

            x[i] = x[i] - 2*Dx

            N = f(x)
            
            #Ako je vrijednost i u drugom smjeru gora nego u početnoj, vrati točku na početnu
            if N>P:

                x[i] = x[i] + Dx

    return x

def HookJeevesTrazenjeUnutarnjeTocke(nejednakosti,x0,vx = 1,e = 10e-6,trace = False):
    Xp = np.array(x0).astype(float)
    Xb = np.array(x0).astype(float)

    while True:


        f = lambda x: unutarnjaTocka(x,nejednakosti)
        #Istrazi najbolju tocku
        Xn = istrazi2(f,Xp.copy(),vx)
        Xn = np.array(Xn).astype(float)

        # Ako je novodobivena tocka bolja od prethodne, prihvati baznu tocku
        if f(Xn) < f(Xb):
            Xp = 2*Xn - Xb
            Xb = Xn
        #Ako novodobivena tocka nije bolja, smanji pomake(Dx)
        else:
            vx = vx/2
            Xp = Xb.copy()

        #Uvjet zaustavljanja
        if vx <= e:
            return Xb


def f1(x):

    return 100*(x[1]-x[0]**2)**2+(1-x[0])**2

f1pox1 = lambda x: 2*(200*x[0]**3-200*x[0]*x[1]+x[0]-1)
f1pox2 = lambda x: 200*(x[1]-x[0]**2)

f1pox1_pox1 = lambda x: -400*(x[1]-x[0]**2)+800*x[0]**2+2
f1pox1_pox2 = lambda x: -400*x[0]
f1pox2_pox1 = lambda x: -400*x[0]
f1pox2_pox2 = lambda x: 200

eksplicitX = lambda x: 1 if x[0] >= -100 and x[0] <= 100 else 0
eksplicitY = lambda x: 1 if x[1] >= -100 and x[1] <= 100 else 0

implicit1 = lambda x: x[1]-x[0]
implicit2 = lambda x: 2-x[0]

def f2(x):

    return (x[0]-4)**2+4*(x[1]-2)**2

f2pox1 = lambda x: 2*(x[0]-4)
f2pox2 = lambda x: 8*(x[1]-2)

f2pox1_pox1 = lambda x: 2
f2pox1_pox2 = lambda x: 0
f2pox2_pox1 = lambda x: 0
f2pox2_pox2 = lambda x: 8

def f3(x):

    return (x[0]-2)**2 + (x[1]+3)**2

f3pox1 = lambda x: 2*(x[0]-2)
f3pox2 = lambda x: 2*(x[1]+3)

def f4(x):

    return (x[0]-3)**2+(x[1])**2

f4nejednakost1 = lambda x: 3-x[0]-x[1]
f4nejednakost2 = lambda x: 3+1.5*x[0]+x[1]
f4jednakost = lambda x: x[1]-1


print("Zadatak 1.\n")

print("Bez računanja optimalnog iznosa koraka:")
functionCalls,gradientCalculations,xmin = gradientDescent(f3,[f3pox1,f3pox2],np.array([0,0]))
print("Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta".format(xmin,functionCalls,gradientCalculations))

print("\nSa računanjem optimalnog iznosa koraka:")
functionCalls,gradientCalculations,xmin = gradientDescent(f3,[f3pox1,f3pox2],np.array([0,0]),precise = True)
print("Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta".format(xmin,functionCalls,gradientCalculations))

print("\nZadatak 2.\n")

print("Funkcija 1:")
functionCalls,gradientCalculations,xmin = gradientDescent(f1,[f1pox1,f1pox2],np.array([-1.9,2]),precise = True)
print("Gradijentnim spustom -> Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta".format(xmin,functionCalls,gradientCalculations))
functionCalls,gradientCalculations,hCalculations,xmin = Newton_Raphson(f1, [f1pox1,f1pox2],[[f1pox1_pox1,f1pox1_pox2],[f1pox2_pox1,f1pox2_pox2]],np.array([-1.9,2]),precise = True)
print("NewtonRaphson -> Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta, H racunan {} puta".format(xmin,functionCalls,gradientCalculations,hCalculations))

print("\nFunkcija 2:")
functionCalls,gradientCalculations,xmin = gradientDescent(f2,[f2pox1,f2pox2],np.array([0.1,0.3]),precise = True)
print("Gradijentnim spustom -> Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta".format(xmin,functionCalls,gradientCalculations))
functionCalls,gradientCalculations,hCalculations,xmin = Newton_Raphson(f2, [f2pox1,f2pox2],[[f2pox1_pox1,f2pox1_pox2],[f2pox2_pox1,f2pox2_pox2]],np.array([0.1,0.3]),precise = True)
print("NewtonRaphson -> Minimum u {}, funkcija pozvana {} puta, gradijent racunan {} puta, H racunan {} puta".format(xmin,functionCalls,gradientCalculations,hCalculations))

print("\nZadatak 3.\n")
minimum = boxovPostupak(f1,[eksplicitX,eksplicitY],[implicit1,implicit2],2,np.array([-1.9,2]))
print("Minimum funkcije 1 uz ograničenja iznosi: {}".format(minimum))
minimum = boxovPostupak(f2,[eksplicitX,eksplicitY],[implicit1,implicit2],2,np.array([0.1,0.3]))
print("Minimum funkcije 2 uz ograničenja iznosi: {}".format(minimum))

print("\nZadatak 4.\n")

minimum = HookJeevesUzOgranicenja(f1,[implicit1,implicit2],[],np.array([-1.9,2]))
print("Minimum funkcije 1 uz ograničenja iznosi: {}".format(minimum))
minimum = HookJeevesUzOgranicenja(f1,[implicit1,implicit2],[],np.array([0.5,2]))
print("Minimum funkcije 1 (uz novu početnu točku [0.5,2]) uz ograničenja iznosi: {}".format(minimum))
minimum = HookJeevesUzOgranicenja(f2,[implicit1,implicit2],[],np.array([0.1,0.3]))
print("Minimum funkcije 2 uz ograničenja iznosi: {}".format(minimum))

print("\nZadatak 5.\n")

minimum = HookJeevesUzOgranicenja(f4,[f4nejednakost1,f4nejednakost2],[f4jednakost],np.array([5,5]))
print("Minimum funkcije 4 uz ograničenja iznosi: {}".format(minimum))
tocka = HookJeevesTrazenjeUnutarnjeTocke([f4nejednakost1,f4nejednakost2],np.array([5,5]))
print("Pronađena unutarnja točka: {}".format(tocka))
minimum = HookJeevesUzOgranicenja(f4,[f4nejednakost1,f4nejednakost2],[f4jednakost],tocka)
print("Minimum funkcije 4 uz ograničenja iznosi: {}".format(minimum))
