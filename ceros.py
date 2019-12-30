import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
#%%
def function(x):
    return np.sin(x)#sp.eval_legendre(6,x)
a = -1.0
b = 7.0
#%%
x = np.linspace(a,b,100)
plt.plot(x,function(x))
plt.title("$Polinomio\ de\ Legendre\ m=6$")
plt.grid(True)
plt.axhline(y=0, c='k')
plt.axvline(x=0, c='k')
plt.show()
#%%

# Encontrar ceros mediante el cambio de signo de la funcion f en el intervalo (a,b) con paso dx
def signo(f,a,b,dx):
    x = a
    sgn = f(x) > 0
    ceros = []
    while x < b:
        if f(x) == 0:
            ceros.append(x) #Encontramos un cero
        else:
            if (f(x) > 0) != sgn:
                ceros.append(x) #Nos pasamos mi rey
        sgn = f(x) > 0
        x += dx
    return ceros

# Encontrar ceros mediante el metodo de biseccion de la funcion f cerca de los puntos en ceros y con tolerancia tol
def biseccion(f,zeros,tol):
    new_zeros = []
    pasos = []
    for cero in zeros:
        #Definimos un nuevo intervalo de ancho 0.1 cerca de cada cero.
        a = cero - 0.05
        b = cero + 0.05
        if f(a)*f(b) > 0:
            raise Exception('Esta vuelta no funciona si no hay cambio de signo')
        c = (a+b)/2.0
        paso = 0
        while (b-a)/2.0 > tol:
            if f(c) == 0:
                new_zeros.append(c)
            else:
                if f(a)*f(c) < 0:
                    b = c
                else:
                    a = c
            c = (a+b)/2.0
            paso += 1
        new_zeros.append(c)
        pasos.append(paso)
    return new_zeros, pasos

def derivative(f, x, h):
    return (f(x+h) - f(x-h))/(2.0*h)

def Newton_Raphson(f, f_prime, h, zeros,tol):
    new_zeros = []
    pasos = []
    for cero in zeros:
        c = cero - f(cero)/f_prime(f,cero,h)
        count = 1
        while abs(f(c)) > tol:
            c = c - f(c)/f_prime(f,c,h)
            count += 1
        new_zeros.append(c)
        pasos.append(count)
    return new_zeros, pasos

#%%
dx = 0.001
ceros_sin = signo(function, a,b, dx)
ceros_sin

ceros_bisec, pasos = biseccion(function,ceros_sin,0.00001)
ceros_bisec
pasos

ceros_newton, pasos_newton = Newton_Raphson(function, derivative, 0.00001, ceros_sin, 0.000001)
ceros_newton
pasos_newton
