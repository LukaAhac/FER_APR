from math import sqrt
from math import sin
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


#print(zlatniRez(lambda x : pow((2+x),2)+4,tocka = 0))

def poKoOsima(f,x0,n,e=10E-6):

    x = np.array(x0)
    emptyVector = [0]*n
    fCallSum = 0

    while True:

        xs = x.copy()

        for i in range(0,n):

            #Sastavi ei vektor
            ei = emptyVector.copy()
            ei[i] = 1
            ei = np.array(ei)

            left,right,fCall = zlatniRez(lambda lamb: f(x+lamb*ei),e=e,tocka=0)
            fCallSum += fCall
            centar = (left+right)/2

            x = x + centar*ei

            if not(np.linalg.norm(x-xs)>e):
                return [list(x),fCallSum]

#rez = poKoOsima(lambda x: x[0]*x[0]+x[1]*x[1], [2,2],2)

#print(rez)

def simpleks(f,x0,e = 10E-6,alfa=1,beta=0.5,gamma=2,sigma=0.5,pomak=1,trace=False):

    tockeSimpleksa = []
    tockeSimpleksa.append(np.array(x0))

    #Racunanje tocaka simpleksa
    for index in range(0,len(x0)):
        tocka = x0.copy()
        tocka[index] += pomak
        tockeSimpleksa.append(np.array(tocka))

        fCall = 0
        korak = 0

    while True:
        
        korak += 1

        h = l = 0
        fXh = fXl = f(tockeSimpleksa[0])

        fCall += 1

        #Određiavnje h i l
        for index in range(1,len(tockeSimpleksa)):

            value = f(tockeSimpleksa[index])

            fCall += 1

            if value > fXh:
                fXh = value
                h = index
            if value < fXl:
                fXl = value
                l = index
                
        Xh = tockeSimpleksa[h]
        Xl = tockeSimpleksa[l]
        #Određivanje centroida

        Xc = [0]*len(x0)
        Xc = np.array(Xc)

        for index in range(0,len(tockeSimpleksa)):
            if index != h:
                Xc = Xc + tockeSimpleksa[index]

        Xc = Xc*1/(len(tockeSimpleksa)-1)

        fXc = f(Xc)
        
        if trace == True:
            print("Korak {} ->Centroid je u:{}, a vrijednost funkcije u njemu je: {}".format(korak,Xc,fXc))

        #Refleksija
        Xr = (1+alfa)*Xc-alfa*Xh
        fXr = f(Xr)
        fCall += 1

        #Ako je vrijednost u reflektiranoj tocki manja od najmanje u centroidu
        if fXr < fXl:
            #Ekspanzija
            Xe = Xc + gamma*(Xr-Xc)
            fXe = f(Xe)
            fCall += 1

            #Ako je vrijednost u ekspandiranoj manja od najmanje u centroidu zamjeni Xh sa njom
            if fXe < fXl:
                Xh = Xe
                tockeSimpleksa[h] = Xe
            #Inace zamjeni Xh sa reflektiranom
            else:
                Xh = Xr
                tockeSimpleksa[h] = Xr
        else:

            #Ispitaj da li je vrijednost u reflektiranoj veca od svih vrijednosti osim u najgoroj Xh
            condition = True
            for index in range(0,len(tockeSimpleksa)):

                if index != h:
                    fCall += 1
                    if not(fXr > f(tockeSimpleksa[index])):
                        condition = False
            #Ako jest
            if condition:
                #Ako je bolja samo od najgore zamjeni najgoru njome
                if fXr < fXh:
                    Xh = Xr
                    fXh = f(Xh)
                    fCall += 1
                    tockeSimpleksa[h] = Xr

                #Kontrakcija
                Xk = (1-beta)*Xc + beta*Xh
                fXk = f(Xk)
                fCall += 1

                #Ako je kontrahirana tocka bolja od najgore zamjeni najgoru njome
                if fXk < fXh:
                    Xh = Xk
                    tockeSimpleksa[h] = Xk
                #Ako nije pomakni sve prema Xl
                else:
                    for index in range(0,len(tockeSimpleksa)):
                        tockeSimpleksa[index] = sigma*(tockeSimpleksa[index]+tockeSimpleksa[l])

            else:
                tockeSimpleksa[h] = Xr

        uvjetZaustavljanja = 0
        for index in range(0,len(tockeSimpleksa)):
            uvjetZaustavljanja += pow((f(tockeSimpleksa[index])-fXc),2)

        uvjetZaustavljanja /= len(tockeSimpleksa)

        if uvjetZaustavljanja<e:

            maxf = f(tockeSimpleksa[0])
            maxi = 0
            for index in range(1,len(tockeSimpleksa)):
                value = f(tockeSimpleksa[index])
                if value<maxf:
                    maxf = value
                    maxi = index

            return [tockeSimpleksa[maxi],fCall]



#print(simpleks(lambda x: (x[0]-2)*(x[0]-2)+x[1]*x[1],[2,0],trace=True))

def istrazi(f, Xp, Dx):

    x = Xp
    fCall = 0

    #Iteriraj po komponentama vektora x
    for i in range(0,len(Xp)):
        P = f(x)
        fCall+= 1

        #Povećaj i-tu komponentu za Dx
        x[i] = x[i] + Dx

        N = f(x)
        fCall+= 1

        #Ako je vrijednost u novodobivenoj točki gora nego u početnoj, idi u drugom smjeru
        if N>P:

            x[i] = x[i] - 2*Dx

            N = f(x)
            fCall+= 1
            
            #Ako je vrijednost i u drugom smjeru gora nego u početnoj, vrati točku na početnu
            if N>P:

                x[i] = x[i] + Dx

    return [x,fCall]

def HookJeeves(f,x0,vx = 1,e = 10e-6,trace = False):
    Xp = np.array(x0).astype(float)
    Xb = np.array(x0).astype(float)
    fCall = 0
    korak = 0

    while True:

        korak += 1
        #Istrazi najbolju tocku
        Xn,additionalFCall = istrazi(f,Xp.copy(),vx)
        Xn = np.array(Xn).astype(float)

        if trace == True:
            print("Korak: {}, Xb: {}, Xp: {}, Xn: {}".format(korak,Xb,Xp,Xn))


        fCall += additionalFCall

        # Ako je novodobivena tocka bolja od prethodne, prihvati baznu tocku
        fCall += 2
        if f(Xn) < f(Xb):
            Xp = 2*Xn - Xb
            Xb = Xn
        #Ako novodobivena tocka nije bolja, smanji pomake(Dx)
        else:
            vx = vx/2
            Xp = Xb.copy()

        #Uvjet zaustavljanja
        if vx <= e:
            return [Xb,fCall]

#print(HookJeeves(lambda x: x[0]*x[0]+4*x[1]*x[1],[7,3],e=0.25,trace = True))

def f1(x):
    return 100*pow(x[1]-pow(x[0],2),2)+pow((1-x[0]),2)

def f2(x):
    return pow((x[0]-4),2)+4*pow((x[1]-2),2)

def f3(x):

    rez = 0
    for i in range(0,len(x)):
        rez += pow(x[i]-3,2)
    return rez

def f3_2(x):
    rez = 0
    for i in range(0,len(x)):
        rez += pow((x[i]-i-1),2)
    return rez

def f4(x):
    return abs((x[0]-x[1])*(x[0]+x[1]))+sqrt(x[0]**2+x[1]**2)

def f6(x):

    suma=0
    
    for i in range(0,len(x)):
        suma += x[i]**2

    return 0.5 + (pow(sin(sqrt(suma)),2)-0.5)/pow((1+0.001*suma),2)



#Zadatak 1
print("Zadatak 1:")

print("")
x0 = 10
print("Početna točka: {}\n".format(x0))

zlatniRez1 = zlatniRez(lambda x : pow((x-3),2),tocka = x0)
print("Zlatni rez - interval: [{},{}], broj poziva funkcije: {}".format(zlatniRez1[0],zlatniRez1[1],zlatniRez1[2]))
poOsima1 = poKoOsima(f3,[x0],1)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f3,[x0])
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f3,[x0])
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))


print("")
x0 = 50
print("Početna točka: {}\n".format(x0))

zlatniRez1 = zlatniRez(lambda x : pow((x-3),2),tocka = x0)
print("Zlatni rez - interval: [{},{}], broj poziva funkcije: {}".format(zlatniRez1[0],zlatniRez1[1],zlatniRez1[2]))
poOsima1 = poKoOsima(f3,[x0],1)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f3,[x0])
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f3,[x0])
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("")
x0 = 100
print("Početna točka: {}\n".format(x0))

zlatniRez1 = zlatniRez(lambda x : pow((x-3),2),tocka = x0)
print("Zlatni rez - interval: [{},{}], broj poziva funkcije: {}".format(zlatniRez1[0],zlatniRez1[1],zlatniRez1[2]))
poOsima1 = poKoOsima(f3,[x0],1)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f3,[x0])
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f3,[x0])
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("\nZadatak 2:\n")

print("Funkcija 1:\n")
f1x0 = [-1.9,2]

poOsima1 = poKoOsima(f1,f1x0,2)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f1,f1x0)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f1,f1x0)
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))


print("\nFunkcija 2:\n")
f2x0 = [0.1,0.3]

poOsima1 = poKoOsima(f2,f2x0,2)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f2,f2x0)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f2,f2x0)
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("\nFunkcija 3:\n")
f3x0 = [0,0,0,0,0]

poOsima1 = poKoOsima(f3_2,f3x0,5)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f3_2,f3x0)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f3_2,f3x0)
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("\nFunkcija 4:\n")
f4x0 = [5.1,1.1]

poOsima1 = poKoOsima(f4,f4x0,2)
print("Traženje po koordinatnim osima - Točka minimuma: {}, broj poziva funkcije: {}".format(poOsima1[0],poOsima1[1]))
simpleks1 = simpleks(f4,f4x0)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f4,f4x0)
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("\nZadatak 3:")

print("\nFunkcija 4:\n")
f4x0 = [5,5]

simpleks1 = simpleks(f4,f4x0)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
HJmetoda1 = HookJeeves(f4,f4x0)
print("HookJeeves metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(list(HJmetoda1[0]),HJmetoda1[1]))

print("\nZadatak 4:")

print("\nFunkcija 1, x0 = [0.5,0.5]:\n")
f1x0 = [0.5,0.5]

print("Korak = 1")
simpleks1 = simpleks(f1,f1x0,pomak = 1)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 5")
simpleks1 = simpleks(f1,f1x0,pomak = 5)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 10")
simpleks1 = simpleks(f1,f1x0,pomak = 10)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 20")
simpleks1 = simpleks(f1,f1x0,pomak = 20)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))

print("\nFunkcija 1, x0 = [20,20]:\n")
f1x0 = [20,20]

print("Korak = 1")
simpleks1 = simpleks(f1,f1x0,pomak = 1)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 5")
simpleks1 = simpleks(f1,f1x0,pomak = 5)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 10")
simpleks1 = simpleks(f1,f1x0,pomak = 10)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))
print("Korak = 20")
simpleks1 = simpleks(f1,f1x0,pomak = 20)
print("Simpleks metoda - Točka minimuma: {}, broj poziva funkcije: {}".format(simpleks1[0],simpleks1[1]))


print("\nZadatak 5:\n")

successfull = 0
numberOfDots = 10000

for i in range(0,numberOfDots):

    tocka = np.random.rand(2)
    tocka = tocka*100
    tocka = tocka-50

    minimum = HookJeeves(f6,x0 = tocka)

    if f6(minimum[0]) < 10e-4:
        successfull += 1

print("Uspješnot pronalaženja globalnog minimuma funkcije f6 HookJeeves metodom iznosi približno {}%".format(successfull/numberOfDots*100))