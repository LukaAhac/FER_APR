import numpy as np
from math import sin,cos
import matplotlib.pyplot as plt

def euler(A,B,rt_fun,t0, T, t_start, t_max):

    values = list()

    current_solution = t0

    while t_start < t_max:
        
        #print(current_solution)
        
        rt = rt_fun(t_start)

        current_solution = current_solution + T * (np.matmul(A,current_solution) + np.matmul(B,rt))

        values.append(current_solution)

        t_start += T

    return values

def reverse_euler(A,B,rt_fun,t0, T, t_start, t_max):

    current_solution = t0
    values = list()

    while t_start < t_max:

        #print(current_solution)

        rt = rt_fun(t_start)

        current_solution = np.matmul(np.linalg.inv(np.identity(A.shape[0]) - A * T),current_solution) + np.matmul(np.matmul(np.linalg.inv(np.identity(A.shape[0]) - A * T)*T,B),rt)

        values.append(current_solution)

        t_start += T

    return values

def trapezni(A,B,rt,t0,T,t_start,t_max):

    current_solution = t0
    values = list()

    while t_start < t_max:

        #print(current_solution)

        rt = rt_fun(t_start)
        rt_next = rt_fun(t_start+T)

        current_solution = np.matmul(np.matmul(np.linalg.inv(np.identity(A.shape[0]) - A * T/2),(np.identity(A.shape[0]) + A * T/2)),current_solution) + np.matmul(np.matmul(np.linalg.inv(np.identity(A.shape[0]) - A * T/2)*T/2,B),rt+rt_next)

        values.append(current_solution)

        t_start += T

    return values

def runge_kutta(A,B,rt,t0,T,t_start,t_max):

    def value(x,t,A,B,rt_fun):

        return np.matmul(A,x) + np.matmul(B,rt_fun(t))

    current_solution = t0
    values = list()

    while t_start < t_max:

        #print(current_solution)
        
        m1 = T * value(current_solution,t_start, A, B, rt_fun)
        m2 = T * value(current_solution + m1/2, t_start + T/2, A, B, rt_fun)
        m3 = T * value(current_solution + m2/2, t_start + T/2, A, B, rt_fun)
        m4 = T * value(current_solution + m3, t_start + T, A, B, rt_fun)

        current_solution = current_solution + 1/6 * (m1 + 2 * m2 + 2 * m3 + m4)

        values.append(current_solution)

        t_start += T

    return values

def pece(A,B,rt,t0,T,t_start,t_max):

    current_solution = t0
    values = list()

    while t_start < t_max:

        rt = rt_fun(t_start)
        rt_next = rt_fun(t_start+T)

        #euler
        x_next =  current_solution + T * (np.matmul(A,current_solution) + np.matmul(B,rt))

        #trapez

        current_solution = current_solution + T/2 * (np.matmul(A,current_solution) + np.matmul(B,rt) + np.matmul(A,x_next) + np.matmul(B,rt_next))

        values.append(current_solution)

        t_start += T

    return values

def pece2(A,B,rt,t0,T,t_start,t_max):

    current_solution = t0
    values = list()

    while t_start < t_max:

        rt = rt_fun(t_start)

        #euler
        x_next = current_solution + T * (np.matmul(A,current_solution) + np.matmul(B,rt))

        #Obrnuti euler
        x_next = current_solution + T * (np.matmul(A,x_next) + np.matmul(B,rt))
        current_solution = current_solution + T * (np.matmul(A,x_next) + np.matmul(B,rt))

        values.append(current_solution)

        t_start += T

    return values


#ZADATAK 1.
print("1. Zadatak\n")

A = np.array([[0,1],[-1,0]])
B = np.array([[0,0],[0,0]])
rt_fun = lambda x: np.array([[0],[0]])
t0 = np.array([[1],[1]])

x1 = lambda x: cos(x) + sin(x)
x2 = lambda x: cos(x) - sin(x)

x1_values = []
x2_values = []

for x in range(1,1001):

    x1_values.append(x1(x/100))
    x2_values.append(x2(x/100))

euler_values = euler(A,B,rt_fun,t0,0.01,0,10)
reverse_euler_values = reverse_euler(A,B,rt_fun,t0,0.01,0,10)
trapez_values = trapezni(A,B,rt_fun,t0,0.01,0,10)
runge_kutta_values = runge_kutta(A,B,rt_fun,t0,0.01,0,10)
pece_values = pece(A,B,rt_fun,t0,0.01,0,10)
pece2_values = pece2(A,B,rt_fun,t0,0.01,0,10)

euler_error = [0,0]
reverse_euler_error = [0,0]
trapez_error = [0,0]
runge_kutta_error = [0,0]
pece_error = [0,0]
pece2_error = [0,0]

for x in range(0,1000):

    euler_error[0] += abs(euler_values[x][0]-x1_values[x])
    euler_error[1] += abs(euler_values[x][1]-x2_values[x])
    reverse_euler_error[0] += abs(reverse_euler_values[x][0]-x1_values[x])
    reverse_euler_error[1] += abs(reverse_euler_values[x][1]-x2_values[x])
    trapez_error[0] += abs(trapez_values[x][0]-x1_values[x])
    trapez_error[1] += abs(trapez_values[x][1]-x2_values[x])
    runge_kutta_error[0] += abs(runge_kutta_values[x][0]-x1_values[x])
    runge_kutta_error[1] += abs(runge_kutta_values[x][1]-x2_values[x])
    pece_error[0] += abs(pece_values[x][0]-x1_values[x])
    pece_error[1] += abs(pece_values[x][1]-x2_values[x])
    pece2_error[0] += abs(pece2_values[x][0]-x1_values[x])
    pece2_error[1] += abs(pece2_values[x][1]-x2_values[x])

print("Pogreška kod Eulera: x1 - {}, x2 - {}".format(euler_error[0],euler_error[1]))
print("Pogreška kod obrnutog Eulera: x1 - {}, x2 - {}".format(reverse_euler_error[0],reverse_euler_error[1]))
print("Pogreška kod trapeza: x1 - {}, x2 - {}".format(trapez_error[0],trapez_error[1]))
print("Pogreška kod Runge-Kutta: x1 - {}, x2 - {}".format(runge_kutta_error[0],runge_kutta_error[1]))
print("Pogreška kod PECE: x1 - {}, x2 - {}".format(pece_error[0],pece_error[1]))
print("Pogreška kod PECE2: x1 - {}, x2 - {}".format(pece2_error[0],pece2_error[1]))

x_axis = list(range(1001))
x_axis = np.array(x_axis)
x_axis = x_axis / 100

x1_euler = []
x2_euler = []
x1_reverse_euler = []
x2_reverse_euler = []
x1_trapez = []
x2_trapez = []
x1_runge_kutta = []
x2_runge_kutta = []
x1_pece = []
x2_pece = []
x1_pece2 = []
x2_pece2 = []

for x in euler_values:

    x1_euler.append(x[0])
    x2_euler.append(x[1])

for x in reverse_euler_values:

    x1_reverse_euler.append(x[0])
    x2_reverse_euler.append(x[1])

for x in trapez_values:

    x1_trapez.append(x[0])
    x2_trapez.append(x[1])

for x in runge_kutta_values:

    x1_runge_kutta.append(x[0])
    x2_runge_kutta.append(x[1])

for x in pece_values:

    x1_pece.append(x[0])
    x2_pece.append(x[1])

for x in pece2_values:

    x1_pece2.append(x[0])
    x2_pece2.append(x[1])

plt.figure("Zadatak1")

plt.subplot(211)
plt.plot(x_axis,x1_euler, label = "Euler")
plt.plot(x_axis,x1_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x1_trapez, label = "Trapez")
plt.plot(x_axis,x1_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x1_pece, label = "PECE")
plt.plot(x_axis,x1_pece2, label = "PECE2")
plt.legend()

plt.subplot(212)
plt.plot(x_axis,x2_euler, label = "Euler")
plt.plot(x_axis,x2_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x2_trapez, label = "Trapez")
plt.plot(x_axis,x2_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x2_pece, label = "PECE")
plt.plot(x_axis,x2_pece2, label = "PECE2")
plt.legend()
plt.show()

#ZADATAK 2.
print("2. Zadatak\n")

A = np.array([[0,1],[-200,-102]])
B = np.array([[0,0],[0,0]])
rt_fun = lambda x: np.array([[0],[0]])
t0 = np.array([[1],[-2]])

#RUNGE KUTA JE STABILAN ZA 0.01 I 0.02, DOK JE ZA 0.03 NESTABILAN

euler_values = euler(A,B,rt_fun,t0,0.01,0,1)
reverse_euler_values = reverse_euler(A,B,rt_fun,t0,0.01,0,1)
trapez_values = trapezni(A,B,rt_fun,t0,0.01,0,1)
runge_kutta_values = runge_kutta(A,B,rt_fun,t0,0.01,0,1)
pece_values = pece(A,B,rt_fun,t0,0.01,0,1)
pece2_values = pece2(A,B,rt_fun,t0,0.01,0,1)

x_axis = list(range(100))
x_axis = np.array(x_axis)
x_axis = x_axis / 100

x1_euler = []
x2_euler = []
x1_reverse_euler = []
x2_reverse_euler = []
x1_trapez = []
x2_trapez = []
x1_runge_kutta = []
x2_runge_kutta = []
x1_pece = []
x2_pece = []
x1_pece2 = []
x2_pece2 = []

for x in euler_values:

    x1_euler.append(x[0])
    x2_euler.append(x[1])

for x in reverse_euler_values:

    x1_reverse_euler.append(x[0])
    x2_reverse_euler.append(x[1])

for x in trapez_values:

    x1_trapez.append(x[0])
    x2_trapez.append(x[1])

for x in runge_kutta_values:

    x1_runge_kutta.append(x[0])
    x2_runge_kutta.append(x[1])

for x in pece_values:

    x1_pece.append(x[0])
    x2_pece.append(x[1])

for x in pece2_values:

    x1_pece2.append(x[0])
    x2_pece2.append(x[1])

plt.figure("Zadatak2")

plt.subplot(211)
plt.plot(x_axis,x1_euler, label = "Euler")
plt.plot(x_axis,x1_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x1_trapez, label = "Trapez")
plt.plot(x_axis,x1_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x1_pece, label = "PECE")
plt.plot(x_axis,x1_pece2, label = "PECE2")
plt.legend()

plt.subplot(212)
plt.plot(x_axis,x2_euler, label = "Euler")
plt.plot(x_axis,x2_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x2_trapez, label = "Trapez")
plt.plot(x_axis,x2_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x2_pece, label = "PECE")
plt.plot(x_axis,x2_pece2, label = "PECE2")
plt.legend()
plt.show()


#ZADATAK 3.
print("3. Zadatak\n")

A = np.array([[0,-2],[1,-3]])
B = np.array([[2,0],[0,3]])
rt_fun = lambda x: np.array([[1],[1]])
t0 = np.array([[1],[3]])

euler_values = euler(A,B,rt_fun,t0,0.01,0,10)
reverse_euler_values = reverse_euler(A,B,rt_fun,t0,0.01,0,10)
trapez_values = trapezni(A,B,rt_fun,t0,0.01,0,10)
runge_kutta_values = runge_kutta(A,B,rt_fun,t0,0.01,0,10)
pece_values = pece(A,B,rt_fun,t0,0.01,0,10)
pece2_values = pece2(A,B,rt_fun,t0,0.01,0,10)

x_axis = list(range(1001))
x_axis = np.array(x_axis)
x_axis = x_axis / 100

x1_euler = []
x2_euler = []
x1_reverse_euler = []
x2_reverse_euler = []
x1_trapez = []
x2_trapez = []
x1_runge_kutta = []
x2_runge_kutta = []
x1_pece = []
x2_pece = []
x1_pece2 = []
x2_pece2 = []

for x in euler_values:

    x1_euler.append(x[0])
    x2_euler.append(x[1])

for x in reverse_euler_values:

    x1_reverse_euler.append(x[0])
    x2_reverse_euler.append(x[1])

for x in trapez_values:

    x1_trapez.append(x[0])
    x2_trapez.append(x[1])

for x in runge_kutta_values:

    x1_runge_kutta.append(x[0])
    x2_runge_kutta.append(x[1])

for x in pece_values:

    x1_pece.append(x[0])
    x2_pece.append(x[1])

for x in pece2_values:

    x1_pece2.append(x[0])
    x2_pece2.append(x[1])

plt.figure("Zadatak3")

plt.subplot(211)
plt.plot(x_axis,x1_euler, label = "Euler")
plt.plot(x_axis,x1_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x1_trapez, label = "Trapez")
plt.plot(x_axis,x1_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x1_pece, label = "PECE")
plt.plot(x_axis,x1_pece2, label = "PECE2")
plt.legend()

plt.subplot(212)
plt.plot(x_axis,x2_euler, label = "Euler")
plt.plot(x_axis,x2_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x2_trapez, label = "Trapez")
plt.plot(x_axis,x2_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x2_pece, label = "PECE")
plt.plot(x_axis,x2_pece2, label = "PECE2")
plt.legend()
plt.show()


#ZADATAK 4.
print("4. Zadatak\n")

A = np.array([[1,-5],[1,-7]])
B = np.array([[5,0],[0,3]])
rt_fun = lambda x: np.array([[x],[x]])
t0 = np.array([[-1],[3]])

euler_values = euler(A,B,rt_fun,t0,0.01,0,1)
reverse_euler_values = reverse_euler(A,B,rt_fun,t0,0.01,0,1)
trapez_values = trapezni(A,B,rt_fun,t0,0.01,0,1)
runge_kutta_values = runge_kutta(A,B,rt_fun,t0,0.01,0,1)
pece_values = pece(A,B,rt_fun,t0,0.01,0,1)
pece2_values = pece2(A,B,rt_fun,t0,0.01,0,1)

x_axis = list(range(100))
x_axis = np.array(x_axis)
x_axis = x_axis / 100

x1_euler = []
x2_euler = []
x1_reverse_euler = []
x2_reverse_euler = []
x1_trapez = []
x2_trapez = []
x1_runge_kutta = []
x2_runge_kutta = []
x1_pece = []
x2_pece = []
x1_pece2 = []
x2_pece2 = []

for x in euler_values:

    x1_euler.append(x[0])
    x2_euler.append(x[1])

for x in reverse_euler_values:

    x1_reverse_euler.append(x[0])
    x2_reverse_euler.append(x[1])

for x in trapez_values:

    x1_trapez.append(x[0])
    x2_trapez.append(x[1])

for x in runge_kutta_values:

    x1_runge_kutta.append(x[0])
    x2_runge_kutta.append(x[1])

for x in pece_values:

    x1_pece.append(x[0])
    x2_pece.append(x[1])

for x in pece2_values:

    x1_pece2.append(x[0])
    x2_pece2.append(x[1])

plt.figure("Zadatak4")

plt.subplot(211)
plt.plot(x_axis,x1_euler, label = "Euler")
plt.plot(x_axis,x1_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x1_trapez, label = "Trapez")
plt.plot(x_axis,x1_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x1_pece, label = "PECE")
plt.plot(x_axis,x1_pece2, label = "PECE2")
plt.legend()

plt.subplot(212)
plt.plot(x_axis,x2_euler, label = "Euler")
plt.plot(x_axis,x2_reverse_euler, label = "Reverse Euler")
plt.plot(x_axis,x2_trapez, label = "Trapez")
plt.plot(x_axis,x2_runge_kutta, label = "Range-Kutta")
plt.plot(x_axis,x2_pece, label = "PECE")
plt.plot(x_axis,x2_pece2, label = "PECE2")
plt.legend()
plt.show()