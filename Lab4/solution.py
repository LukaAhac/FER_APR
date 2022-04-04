from math import ceil
from math import floor
from math import log
from math import sin
from math import sqrt
from random import randint
from random import uniform
from random import sample
from random import choice
from statistics import median


class f1():

    def __init__(self):
        self.calls = 0

    def calc(self,x):

        self.calls += 1

        return 100*(x[1]-x[0]**2)**2 + (1-x[0])**2


class f2():

    def __init__(self):
        self.calls = 0

    def calc(self,x):

        self.calls += 1

        result = 0

        for index in range(len(x)):
            result += (x[index]-(index+1))**2
        
        return result

class f3():

    def __init__(self):
        self.calls = 0

    def calc(self,x):

        self.calls += 1

        result = 0

        for index in range(len(x)):
            result += x[index]**2
        
        return 0.5 + ((sin(sqrt(result)))**2-0.5)/(1+0.001*result)**2

class f4():

    def __init__(self):
        self.calls = 0

    def calc(self,x):

        self.calls += 1

        result = 0

        for index in range(len(x)):
            result += x[index]**2
        
        return (result**0.25)*(1+(sin(50*(result**0.1)))**2)

        


class k_turnir_elimination_GA():

    def __init__(self, k = 3, display = "bin", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 10, pm = 0.1, evaluations = 5000, function = None, crossoverType = 1):

        self.k = k
        self.display = display
        self.dg = dg
        self.gg = gg
        self.variables = variables
        self.precision = precision
        self.popSize = popSize
        self.pm = pm
        self.evaluations = evaluations
        self.function = function
        self.crossoverType = crossoverType

        if self.display == "bin":

            self.numberOfBits = self.variables * ceil((log(floor((1+(self.gg-self.dg)*pow(10,self.precision))))/(log(2))))

    class Chromosome:

        def __init__(self, display, variables, dg, gg, size = None):

            self.display = display
            self.size = size
            self.variables = variables
            self.dg = dg
            self.gg = gg

            if display == "bin":
                self.value = []

                for _ in range(size):
                    self.value.append(randint(0,1))

            elif display == "dec":
                self.value = []

                for _ in range(variables):
                    self.value.append(uniform(dg,gg))
            else:
                raise Exception("Wrong display!")

        def decode(self):

            decoded_vars = []
            to_Decode = []
            binaries = []

            for index in range(self.variables):
                to_Decode.append(self.value[int(index*(self.size/self.variables)):int((index+1)*(self.size/self.variables))])

            for binary in to_Decode:
                b = 1
                value = 0
                binary.reverse()
                for bit in binary:
                    value += b*bit
                    b *= 2
                binaries.append(value)

            for b in binaries:
                decoded_vars.append(self.dg+ b/(pow(2,self.size/self.variables)-1)*(self.gg-self.dg))

            return decoded_vars

        def evaluate(self,func):

            if self.display == "bin":
                self.fitnessValue = func.calc(self.decode())
            elif self.display == "dec":
                self.fitnessValue = func.calc(self.value)
            else:
                raise Exception("Wrong display!")

            


    def crossover(self,parent1,parent2,child):

        if self.display == "bin":
            #Binarno s jednom točkom prekida
            if self.crossoverType == 1:
                value = []
                stop = randint(1,self.numberOfBits)
                for index in range(self.numberOfBits):
                    if index < stop:
                        value.append(parent1.value[index])
                    else:
                        value.append(parent2.value[index])
                child.value = value
            #Binarno uniformno
            elif self.crossoverType == 2:
                value = []
                for index in range(self.numberOfBits):
                    value.append(choice([parent1.value[index],parent2.value[index]]))
                child.value = value
            else:
                raise Exception("Wrong crossover type!")
        elif self.display == "dec":
            #Decimalno aritmetičko
            if self.crossoverType == 1:
                value = []
                a = uniform(0,1)
                for index in range(self.variables):
                    value.append(a*parent1.value[index]+(1-a)*parent2.value[index])
                child.value = value
            #Decimalno heurističko
            elif self.crossoverType == 2:
                value = []
                a = uniform(0,1)
                if parent1.fitnessValue > parent2.fitnessValue:
                    for index in range(self.variables):
                        value.append(a*(parent2.value[index]-parent1.value[index])+parent2.value[index])
                else:
                    for index in range(self.variables):
                        value.append(a*(parent1.value[index]-parent2.value[index])+parent1.value[index])
                for index in range(self.variables):
                    if value[index] < self.dg:
                        value[index] = self.dg
                    elif value[index] > self.gg:
                        value[index] = self.gg
                child.value = value

    def mutation(self, unit):

        mutation = uniform(0,1)

        if mutation <= self.pm:

            #Jednolika binarna mutacija
            if self.display == "bin":
                for index in range(self.numberOfBits):
                    if unit.value[index] == 0:
                        unit.value[index] = 1
                    elif unit.value[index] == 1:
                        unit.value[index] = 0
                    else:
                        raise Exception("Illegal value!")

            #Jednolika aritmetička mutacija
            elif self.display == "dec":
                for index in range(self.variables):
                    unit.value[index] = uniform(self.dg,self.gg)

        else:
            pass




    def run(self):

        #stvori populaciju
        self.population = []
        if self.display == "bin":
            for _ in range(self.popSize):
                self.population.append(self.Chromosome(display = self.display, size = self.numberOfBits, variables = self.variables, dg = self.dg, gg = self.gg))
        elif self.display == "dec":
            for _ in range(self.popSize):
                self.population.append(self.Chromosome(display = self.display, variables = self.variables, dg = self.dg, gg = self.gg))
        else:
            raise Exception("Wrong display!")

        #evaluairaj populaciju
        for example in self.population:
            example.evaluate(self.function)

        #ponavljaj dok nije zadovoljen uvjet zaustavljanja
        while(self.function.calls <= self.evaluations):

            #odaberi slucajno k jedniki
            k_samples = sample(self.population,self.k)

            #pronađi najlosiju
            k_samples.sort(key = lambda x: x.fitnessValue, reverse = True)

            #nova = krizanje neke preostale 2
            toCrossover = sample(k_samples[1:],2)
            self.crossover(toCrossover[0],toCrossover[1],k_samples[0])

            #mutacija
            self.mutation(k_samples[0])

            #evaluacija nove
            k_samples[0].evaluate(self.function)
             #dodaj novu u populaciju -> nova samo zamjeni vrijednost stare


        self.population.sort(key = lambda x: x.fitnessValue)

        if self.display == "bin":
            return [self.population[0].decode(),self.population[0].fitnessValue]
        else:
            return [self.population[0].value,self.population[0].fitnessValue]
        


zadatak = 0

if zadatak == 1 or zadatak == 0:
    print("1.Zadatak\n")

    f1dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f1(), crossoverType = 2)
    x,fx = f1dec.run()
    print("Funkcija f1 - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    f1bin = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f1(), crossoverType = 2)
    x,fx = f1bin.run()
    print("Funkcija f1 - BINARNI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    print("\n")

    f2dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = 5, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f2(), crossoverType = 2)
    x,fx = f2dec.run()
    print("Funkcija f2 - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    f2bin = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = 5, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f2(), crossoverType = 1)
    x,fx = f2bin.run()
    print("Funkcija f2 - BINARNI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    print("\n")    

    f3dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f3(), crossoverType = 2)
    x,fx = f3dec.run()
    print("Funkcija f3 - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    f3bin = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f3(), crossoverType = 2)
    x,fx = f3bin.run()
    print("Funkcija f3 - BINARNI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    print("\n")

    f4dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f4(), crossoverType = 2)
    x,fx = f4dec.run()
    print("Funkcija f4 - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    f4bin = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = 2, precision = 3, popSize = 200, pm = 0.3, evaluations = 1000000, function = f4(), crossoverType = 1)
    x,fx = f4bin.run()
    print("Funkcija f4 - BINARNI PRIKAZ: x = {}, f(x) = {}".format(x,fx))
    print("\n")

if zadatak == 2 or zadatak == 0:
    print("Zadatak 2.\n")

    dimensions = [1,3,6,10]

    for index in range(len(dimensions)):

        f3dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = dimensions[index], precision = 3, popSize = 200, pm = 0.3, evaluations = 100000, function = f3(), crossoverType = 2)
        x,fx = f3dec.run()
        print("Funkcija f3 sa {} varijabla - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(dimensions[index],x,fx))

    print("\n")


    for index in range(len(dimensions)):

        f4dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = dimensions[index], precision = 3, popSize = 200, pm = 0.3, evaluations = 100000, function = f4(), crossoverType = 2)
        x,fx = f4dec.run()
        print("Funkcija f3 sa {} varijabla - DECIMLANI PRIKAZ: x = {}, f(x) = {}".format(dimensions[index],x,fx))

    print("\n")

if zadatak == 3 or zadatak == 0:
    print("Zadatak 3\n")

    dimensions = [3,6]

    f3dec3var = []
    f3bin3var = []
    f3dec6var = []
    f3bin6var = []

    f4dec3var = []
    f4bin3var = []
    f4dec6var = []
    f4bin6var = []

    for index in range(len(dimensions)):

        for _ in range(10):

            f3dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = dimensions[index], precision = 4, popSize = 200, pm = 0.3, evaluations = 100000, function = f3(), crossoverType = 2)
            x,fx = f3dec.run()
            if index == 0:
                f3dec3var.append(fx)
            else:
                f3dec6var.append(fx)
            f3dec = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = dimensions[index], precision = 4, popSize = 200, pm = 0.3, evaluations = 100000, function = f3(), crossoverType = 2)
            x,fx = f3dec.run()
            if index == 0:
                f3bin3var.append(fx)
            else:
                f3bin6var.append(fx)

            f4dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = dimensions[index], precision = 4, popSize = 200, pm = 0.3, evaluations = 100000, function = f4(), crossoverType = 2)
            x,fx = f4dec.run()
            if index == 0:
                f4dec3var.append(fx)
            else:
                f4dec6var.append(fx)           
            f4dec = k_turnir_elimination_GA(k = 3, display = "bin", dg = -50, gg = 150, variables = dimensions[index], precision = 4, popSize = 200, pm = 0.3, evaluations = 100000, function = f4(), crossoverType = 2)
            x,fx = f4dec.run()
            if index == 0:
                f4bin3var.append(fx)
            else:
                f4bin6var.append(fx)

    print("Pronađeni minimumi za f3 - DECIMALNO - 3 VAR")
    print(f3dec3var,"\n")
    print("Pronađeni minimumi za f3 - BINARNO - 3 VAR")
    print(f3bin3var,"\n")
    print("Pronađeni minimumi za f3 - DECIMALNO - 6 VAR")
    print(f3dec6var,"\n")
    print("Pronađeni minimumi za f3 - BINARNO - 6 VAR")
    print(f3bin6var,"\n")
    print("Pronađeni minimumi za f4 - DECIMALNO - 3 VAR")
    print(f4dec3var,"\n")
    print("Pronađeni minimumi za f4 - BINARNO - 3 VAR")
    print(f4bin3var,"\n")
    print("Pronađeni minimumi za f4 - DECIMALNO - 6 VAR")
    print(f4dec6var,"\n")
    print("Pronađeni minimumi za f4 - BINARNO - 6 VAR")
    print(f4bin6var,"\n")

if zadatak == 4 or zadatak == 0:
    print("Zadatak4\n")

    popSizes = [30,50,100,200]
    pm = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]

    meidanDict = dict()
    valuesDict = dict()

    for index1 in range(len(popSizes)):

        for index2 in range(len(pm)):

            meidanDict[str(popSizes[index1])+"---"+str(pm[index2])] = 0
            valuesDict[str(popSizes[index1])+"---"+str(pm[index2])] = []


    for index1 in range(len(popSizes)):

        for index2 in range(len(pm)):

            for _ in range(10):


                f3dec = k_turnir_elimination_GA(k = 3, display = "dec", dg = -50, gg = 150, variables = 10, precision = 4, popSize = popSizes[index1], pm = pm[index2], evaluations = 100000, function = f3(), crossoverType = 2)
                x,fx = f3dec.run()

                valuesDict[str(popSizes[index1])+"---"+str(pm[index2])].append(fx)

            meidanDict[str(popSizes[index1])+"---"+str(pm[index2])] = median(valuesDict[str(popSizes[index1])+"---"+str(pm[index2])])

    print(meidanDict)
    print(valuesDict)
    print("\n")

if zadatak == 5 or zadatak == 0:
    print("Zadatak 5\n")

    k_size = [3,5,10,20]

    meidanDict = dict()
    valuesDict = dict()

    for index in range(len(k_size)):

        meidanDict[str(k_size[index])] = 0
        valuesDict[str(k_size[index])] = []

    for index in range(len(k_size)):

        for _ in range(10):

            f3dec = k_turnir_elimination_GA(k = k_size[index], display = "dec", dg = -50, gg = 150, variables = 10, precision = 4, popSize = 200, pm = 0.6, evaluations = 100000, function = f3(), crossoverType = 2)
            x,fx = f3dec.run()   
      
            valuesDict[str(k_size[index])].append(fx)

        meidanDict[str(k_size[index])] = median(valuesDict[str(k_size[index])])


    print(meidanDict)
    print(valuesDict)
    print("\n")
