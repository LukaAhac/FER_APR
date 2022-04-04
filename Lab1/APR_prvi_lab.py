

class Matrix():

    epsilon = 10E-12

    #Constructor that takes path to file where matrix is writen and loads it
    def __init__(self, path=None):

        if path == None:
            self.values = []
            self.rows = 0
            self.columns = 0
            return
        else:
            file = open(path,"r")
            lines = file.read()
            lines = lines.split("\n")

            self.values = []
            self.rows = len(lines)


            for row in lines:
                row = row.split()
                row = list(map(float,row))
                self.values.append(row)

            self.columns = len(row)

    #Method to print on display
    def print(self):

        print("Matrica {}x{}: ".format(self.rows,self.columns))
        for x in range(0,self.rows):
            for y in range(0,self.columns):
                print("{} ".format(self.values[x][y]), end="")
            print("")
        print("")

    #Method to save a matrix into a file, requires path to a file
    def writeToFile(self,path):

        file = open(path,"w")
        for x in range(0,self.rows):
            for y in range(0,self.columns):
                file.write("{} ".format(self.values[x][y]))
            if x != self.rows-1:
                file.write("\n")

    #Addition
    def __add__(self, other):
        
        if self.rows != other.rows or self.columns != other.columns:
            print("Cannot add those!")
            return None

        newMatrix = Matrix()
        newMatrix.rows = self.rows
        newMatrix.columns = self.columns

        for x in range(0,self.rows):
            row = []
            for y in range(0,self.columns):
                row.append(0)
            newMatrix.values.append(row)

        for x in range(0,self.rows):
            for y in range(0,self.columns):
                newMatrix.values[x][y] = self.values[x][y]+other.values[x][y]

        return newMatrix

    #Subtraction
    def __sub__(self, other):
        
        if self.rows != other.rows or self.columns != other.columns:
            print("Cannot sub those!")
            return None

        newMatrix = Matrix()
        newMatrix.rows = self.rows
        newMatrix.columns = self.columns

        for x in range(0,self.rows):
            row = []
            for y in range(0,self.columns):
                row.append(0)
            newMatrix.values.append(row)

        for x in range(0,self.rows):
            for y in range(0,self.columns):
                newMatrix.values[x][y] = self.values[x][y]-other.values[x][y]

        return newMatrix

    #Multiplication
    def __mul__(self, other):
        
        if self.columns != other.rows:
            print("Cannot multiply those!")
            return None

        newMatrix = Matrix()
        newMatrix.rows = self.rows
        newMatrix.columns = other.columns

        for x in range(0,self.rows):
            row = []
            for y in range(0,other.columns):
                row.append(0)
            newMatrix.values.append(row)

        for x in range(0,self.rows):
            for y in range(0,other.columns):
                for k in range(0,self.columns):
                    newMatrix.values[x][y] += self.values[x][k]*other.values[k][y]

        return newMatrix

    #Transpose a matrix
    def __invert__(self):

        newMatrix = Matrix()
        newMatrix.rows = self.columns
        newMatrix.columns = self.rows

        for x in range(0,newMatrix.rows):
            row = []
            for y in range(0,newMatrix.columns):
                row.append(0)
            newMatrix.values.append(row)

        for x in range(0,newMatrix.rows):
            for y in range(0,newMatrix.columns):
                newMatrix.values[x][y] = self.values[y][x]

        return newMatrix

    # method for +=
    def __iadd__(self,other):

        if self.rows != other.rows or self.columns != other.columns:
            print("Cannot add those!")
            return None

        for x in range(0,self.rows):
            for y in range(0,self.columns):
                self.values[x][y] = self.values[x][y]+other.values[x][y]
        
        return self

    # method for -=
    def __isub__(self,other):

        if self.rows != other.rows or self.columns != other.columns:
            print("Cannot sub those!")
            return None

        for x in range(0,self.rows):
            for y in range(0,self.columns):
                self.values[x][y] = self.values[x][y]-other.values[x][y]
        
        return self

    #method for scalar multiplication, performs it on self
    def scalarMul(self,scalar):

        for x in range(0,self.rows):
            for y in range(0,self.columns):
                self.values[x][y] = self.values[x][y]*scalar
        
        return self

    # method for ==
    def __eq__(self,other):

        if self.rows != other.rows or self.columns != other.columns:
            return False
        for x in range(0,self.rows):
            for y in range(0,self.columns):
                if self.values[x][y] != other.values[x][y]:
                    return False
        return True

    #Method that checks if number is close enough to be counted as zero, uses epsiilon defined in Matrix class (Matrix.epsilon)
    def checkZero(self,number):
        if abs(number) < self.epsilon:
            return True
        return False

    #Methods that performs LU decomposition on itself
    def luDecomposition(self):

        if self.rows != self.columns:
            print("Unable to make a LU decomposition to a {}x{} matrux".format(self.rows,self.columns))
            return None

        for i in range(0,self.rows-1):
            for j in range(i+1,self.columns):
                if(self.checkZero(self.values[i][i]) is True):
                    print("Cannot be done - zero division")
                    return None

                self.values[j][i] /= self.values[i][i]
                for k in range(i+1,self.columns):
                    self.values[j][k] -= self.values[j][i]*self.values[i][k]


    #Methods that performs LUP decomposition on itself and it returns P matrix
    def lupDecomposition(self):


        pMatrix = Matrix()
        pMatrix.rows = self.rows
        pMatrix.columns = self.columns
        pMatrix.numberOfSwaps = 0

        for x in range(0,pMatrix.rows):
            row = []
            for y in range(0,pMatrix.columns):
                if x==y:
                    row.append(1)
                else:
                    row.append(0)
            pMatrix.values.append(row)

        if self.rows != self.columns:
            print("Unable to make a LUP decomposition to a {}x{} matrux".format(self.rows,self.columns))
            return None

        for i in range(0,self.rows-1):

            maxRow = i
            for seeMax in range(i,self.rows):
                if abs(self.values[seeMax][i]) > abs(self.values[maxRow][i]):
                    maxRow = seeMax

            if self.checkZero(self.values[maxRow][i]) == True:
                print("Cannot be done - zero division")
                return None

            if maxRow != i:

                pMatrix.numberOfSwaps += 1

                swap = self.values[i]
                self.values[i] = self.values[maxRow]
                self.values[maxRow] = swap

                swap = pMatrix.values[i]
                pMatrix.values[i] = pMatrix.values[maxRow]
                pMatrix.values[maxRow] = swap

            for j in range(i+1,self.columns):
                self.values[j][i] /= self.values[i][i]
                for k in range(i+1,self.columns):
                    self.values[j][k] -= self.values[j][i]*self.values[i][k]

        return pMatrix

    #Method that performs forward supstitution, requres b matrix to be given at least, and p matrix as optional, method returns result matrix (offten later refferd as y matrix)
    def forwardSupstitution(self,bMatrix,pMatrix = None):

        if pMatrix is not None:
            bMatrix = pMatrix*bMatrix

        for i in range(0,self.rows-1):
            for j in range(i+1,self.rows):
                bMatrix.values[j][0] -= self.values[j][i]*bMatrix.values[i][0]

        return bMatrix


    #Method that performs backward supstitution, requres y matrix to be given, method returs result matrix (offten refferd as x matrix)
    def backwardSupstitution(self,yMatrix):

        for i in range(self.rows-1,-1,-1):
            if self.checkZero(self.values[i][i]) == True:
                print("Cannot be done - zero division")
                return None
            yMatrix.values[i][0] /= self.values[i][i]
            for j in range(0,i):
                yMatrix.values[j][0] -= self.values[j][i]*yMatrix.values[i][0]

        return yMatrix

    #method that calcultes and returns an invers
    def inverse(self):

        pMatrix = self.lupDecomposition()

        inverseMatrix = Matrix()
        inverseMatrix.rows = self.rows
        inverseMatrix.columns = self.columns

        for x in range(0,inverseMatrix.rows):
            row = []
            for y in range(0,inverseMatrix.columns):
                row.append(0)
            inverseMatrix.values.append(row)


        for column in range(0,pMatrix.columns):

            e = Matrix()
            e.rows = pMatrix.rows
            e.columns = 1

            for rows in range(0,pMatrix.rows):
                row = []
                row.append(pMatrix.values[rows][column])
                e.values.append(row)


            y = self.forwardSupstitution(e)

            x = self.backwardSupstitution(y)
            
            if x is None:
                return None

            for row in range(0,x.rows):
                inverseMatrix.values[row][column] = x.values[row][0]

        return inverseMatrix

    # method that calculates and returns determinant of a matrix
    def determinant(self):
            
        p = self.lupDecomposition()

        det = pow(-1,p.numberOfSwaps)

        for x in range(0,self.rows):
            det *= self.values[x][x]

        return det



#ZADATCI SA LABOSA

zadatak = 0

if zadatak == 2 or zadatak == 0:
    print("Drugi zadatak\n")

    A2 = Matrix(r"Lab1\Matrices\2A.txt")
    b2 = Matrix(r"Lab1\Matrices\2B.txt")

    print("LU dekompozicijom:\n")
    A2.print()
    print("Nakon dekompozicije:\n")
    A2.luDecomposition()

    print("\n\n")

    A2 = Matrix(r"Lab1\Matrices\2A.txt")
    b2 = Matrix(r"Lab1\Matrices\2B.txt")

    print("LUP dekompozicijom:\n")
    A2.print()
    print("Nakon dekompozicije:\n")
    p = A2.lupDecomposition()
    A2.print()
    y = A2.forwardSupstitution(b2,p)
    print("Vektor y:")
    y.print()
    x = A2.backwardSupstitution(y)
    print("Vektor x:")
    x.print()

if zadatak == 3 or zadatak == 0:
    print("Treći zadatak")

    A3 = Matrix(r"Lab1\Matrices\3A.txt")
    b2 = Matrix(r"Lab1\Matrices\2B.txt")

    print("LU dekompozicijom:\n")
    A3.print()
    print("Nakon dekompozicije:\n")
    A3.luDecomposition()
    #Element (2,2) je 0 pa nema rjesenja
    A3.print()

    y = A3.forwardSupstitution(b2)

    print("Vektor y:")
    y.print()

    print("Vektor x:")
    x = A3.backwardSupstitution(y) #Ne ide radi nule

    print("\n\n")

    A3 = Matrix(r"Lab1\Matrices\3A.txt")
    b2 = Matrix(r"Lab1\Matrices\2B.txt")

    print("LUP dekompozicijom:\n")
    A3.print()
    print("Nakon dekompozicije:\n")
    p = A3.lupDecomposition()
    A3.print()
    y = A3.forwardSupstitution(b2,p)
    print("Vektor y:")
    y.print()
    print("Vektor x:")
    x = A3.backwardSupstitution(y)

if zadatak == 4 or zadatak == 0:
    print("Četvrti zadatak")

    A4 = Matrix(r"Lab1\Matrices\4A.txt")
    b4 = Matrix(r"Lab1\Matrices\4B.txt")

    print("LU dekompozicijom:\n")
    A4.print()
    print("Nakon dekompozicije:\n")
    A4.luDecomposition()
    A4.print()
    y = A4.forwardSupstitution(b4)
    print("Vektor y:")
    y.print()
    print("Vektor x:")
    x = A4.backwardSupstitution(y)
    x.print()
   
    print("\n\n")

    A4 = Matrix(r"Lab1\Matrices\4A.txt")
    b4 = Matrix(r"Lab1\Matrices\4B.txt")

    print("LUP dekompozicijom:\n")
    A4.print()
    print("Nakon dekompozicije:\n")
    p = A4.lupDecomposition()
    A4.print()
    y = A4.forwardSupstitution(b4,p)
    print("Vektor y:")
    y.print()
    print("Vektor x:")
    x = A4.backwardSupstitution(y)
    x.print()

if zadatak == 5 or zadatak == 0:
    print("Peti zadatak")

    A5 = Matrix(r"Lab1\Matrices\5A.txt")
    b5 = Matrix(r"Lab1\Matrices\5B.txt")

    print("LUP dekompozicijom:\n")
    A5.print()
    print("Nakon dekompozicije:\n")
    p = A5.lupDecomposition()
    A5.print()
    y = A5.forwardSupstitution(b5,p)
    print("Vektor y:")
    y.print()
    print("Vektor x:")
    x = A5.backwardSupstitution(y)
    x.print()

if zadatak == 6 or zadatak == 0:
    print("Šesti zadatak")
    
    A6 = Matrix(r"Lab1\Matrices\6A.txt")
    b6 = Matrix(r"Lab1\Matrices\6B.txt")

    print("LUP dekompozicijom:\n")
    A6.print()
    print("Nakon dekompozicije:\n")
    p = A6.lupDecomposition()
    A6.print()
    y = A6.forwardSupstitution(b6,p)
    print("Vektor y:")
    y.print()
    print("Vektor x:")
    x = A6.backwardSupstitution(y)
    x.print()

if zadatak == 7 or zadatak == 0:
    print("Sedmi zadatak")

    A7 = Matrix(r"Lab1\Matrices\7A.txt")
    A7.print()

    print("Inverz: ")
    inverseA7 = A7.inverse()

if zadatak == 8 or zadatak == 0:
    print("Osmi zadatak")

    A8 = Matrix(r"Lab1\Matrices\8A.txt")
    A8.print()

    print("Inverz: ")
    inverseA8 = A8.inverse()
    inverseA8.print()

if zadatak == 9 or zadatak == 0:
    print("Deveti zadatak")

    A9 = Matrix(r"Lab1\Matrices\9A.txt")
    A9.print()

    detA9= A9.determinant()
    print("Determinanta: {}".format(detA9))

if zadatak == 10 or zadatak == 0:
    print("Deseti zadatak")

    A10 = Matrix(r"Lab1\Matrices\10A.txt")
    A10.print()

    detA10= A10.determinant()
    print("Determinanta: {}".format(detA10))