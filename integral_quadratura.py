import numpy as np
import matplotlib.pyplot as plt

#Integra de-inf a inf
#Para integrales con intervalos diferentes, consultar el libro de landau
def funcion(x):
    return np.cos(x)*np.exp(-x*x-x)

def gauss(Npoints,a):
    x = np.zeros((2001))
    w = np.zeros((2001))

    m = i = j = t = t1 = pp = p1 = p2 = p3 = 0

    eps = 3.0E-14

    m = int((Npoints + 1)/2)

    for i in range(1,m+1):
        t = np.cos((np.pi*(i-0.25))/(Npoints + 0.5))
        t1 = 1
        while(abs(t-t1) >=  eps):
            p1 = 1.0
            p2 = 0.0
            for j in range(1,Npoints+1):
                p3 = p2
                p2 = p1
                p1 = ((2.0*j-1) * t * p2 - (j-1) * p3)/j
            pp = Npoints * (t * p1 -p2)/(t*t - 1.0)
            t1 = t
            t = t1 - p1/pp
        x[i-1] = -t
        x[Npoints -i] = t
        w[i-1] = 2.0/((1.0 - t*t)*pp*pp)
        w[Npoints -i] = w[i-1]
    for i in range(0,Npoints):
        w[i] = w[i]*a*(1.0+x[i]*x[i])/(1.0-x[i]*x[i])**2
        x[i] = x[i]*a/(1.0-x[i]*x[i])
    return x, w
def GaussInt(a, eNes, fun):
    integral = []
    for ene in eNes:
        quadra = 0
        x, w = gauss(ene+1, a)
        for i in range(0,ene+1):
            quadra += fun(x[i])*w[i]
        integral.append(quadra)
    return np.array(integral)

N = [20]

Inte_Gauss = GaussInt(1.0, N, funcion)

print("La integral con el metodo de cuadtratura da: " , Inte_Gauss[0])
