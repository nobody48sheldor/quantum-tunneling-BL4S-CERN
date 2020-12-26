import matplotlib.pyplot as plt
import numpy as np
from math import *

a = int(input("render level "))
b = int(input("render level animation = "))

l = float(input("lenght = "))*(10**(-10))
q = l/2.07
r = float(input("size of potential barrer = "))*(10**(-10))
n = int(input("n = "))
h = 6.62607015*(10**(-34))
h_bar = h/(2*np.pi)
m = 9.1093837015*(10**(-31))
E = ((n**2)*(h**2))/(8*m*l**2)
print("the total energy of the particle is ", E)
V_o = float(input("V_o = "))*(1.6*(10**(-19)))
alpha = (np.sqrt(2*m*(V_o-E)))/h_bar
print("alpha ", alpha)

tmax = float(input("tmax = "))*(10**(-6))

xA = np.linspace(0, q, a)
xB = np.linspace(q, (q+r), a)
xC = np.linspace((q+r), l, a)

t = np.linspace(0, tmax, a)


def delta(x):
    s = (np.exp(-alpha*x))/(2*(np.sin(((n*np.pi)/l)*x)))
    return(s)


B = np.sqrt(1/((2*((delta(q))**2)*(q - (l/(2*n*np.pi))*np.sin(((2*n*np.pi)/l)*q))) - (1/(2*alpha))*((np.exp(-alpha*(q+r)))-(np.exp(-alpha*q))) + (2*((delta(q+r))**2)*(l - q + r + (l/(2*n*np.pi))*np.sin(((2*n*np.pi)/l)*(q+r))))))
A = 2*B*delta(q)
C = 2*B*delta(q+r)

I = B**2*(((2*((delta(q))**2)*(q - (l/(2*n*np.pi))*np.sin(((2*n*np.pi)/l)*q))) - (1/(2*alpha))*((np.exp(-alpha*(q+r)))-(np.exp(-alpha*q))) + (2*((delta(q+r))**2)*(l - q + r + (l/(2*n*np.pi))*np.sin(((2*n*np.pi)/l)*(q+r))))))
print("total integral = ", I)

I_ = (B**2*(2*((delta(q+r))**2)*(l - q + r + (l/(2*n*np.pi))*np.sin(((2*n*np.pi)/l)*(q+r)))))*100
print("probability of the particle to be behind the potential barrer = ", I_, "%")


def AP(x):
    s = (A*np.sin(((n*np.pi)/l)*x))**2
    return(s)

def BP(x):
    s = B**2*(np.exp((-2*alpha*x)))
    return(s)

def CP(x):
    s = (C*np.sin(((n*np.pi)/l)*x))**2
    return(s)


plt.plot(xA, AP(xA), color = 'red', label = "density probablity")
plt.plot(xB, BP(xB), color = 'green', label = "density probability in the potential barrer")
plt.plot(xC, CP(xC), color = 'red')

plt.legend()
plt.show()


def AR(x, t):
    s = np.cos(-(E/h_bar)*t) * (A*np.sin(((n*np.pi)/l)*x))
    return(s)

def AC(x, t):
    s = np.sin(-(E/h_bar)*t) * (A*np.sin(((n*np.pi)/l)*x))
    return(s)

def BR(x, t):
    s = np.cos(-(E/h_bar)*t) * (B*(np.exp((-alpha*x))))
    return(s)

def BC(x, t):
    s = np.sin(-(E/h_bar)*t) * (B*(np.exp((-alpha*x))))
    return(s)

def CR(x, t):
    s = np.cos(-(E/h_bar)*t) * (C*np.sin(((n*np.pi)/l)*x))
    return(s)

def CC(x, t):
    s = np.sin(-(E/h_bar)*t) * (C*np.sin(((n*np.pi)/l)*x))
    return(s)

plt.figure()
plt.ion()
plt.plot([], [])

i = 0

while i<a:
    ar = AR(xA, t[i])
    ac = AC(xA, t[i])

    br = BR(xB, t[i])
    bc = BC(xB, t[i])

    cr = CR(xC, t[i])
    cc = CC(xC, t[i])

    print(t[i])

    plt.clf()
    plt.plot(xA, ar, color = 'blue')
    plt.plot(xA, ac, color = 'green')

    plt.plot(xB, br, color = 'blue')
    plt.plot(xB, bc, color = 'green')

    plt.plot(xC, cr, color = 'blue')
    plt.plot(xC, cc, color = 'green')
    plt.xlim([0, l])
    plt.ylim([-(A + A/100), (A + A/100)])
    plt.pause((1/(a*tmax*10**4)))
    i = i + 1
