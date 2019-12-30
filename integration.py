import numpy as np
import matplotlib.pyplot as plt

#%%

#Trapezio
def Decay(x):
    return np.exp(-x)

def PesoT(i, h, n):
    if (i == 1 or i==n):
        w = h/2.0
    else:
        w = h
    return w

def Trapezio(a,b, eNes, fun):
    integral= []
    for ene in eNes:
        suma = 0.0
        h = (b-a)/(ene-1)
        for i in range (1,int(ene)+1):
            suma += PesoT(i,h,ene)*fun(a + (i-1)*h)
        integral.append(suma)
    return np.array(integral)


eNes = [2,10, 20, 40, 80, 160, 320, 640, 1280]


exact = 1.0 - np.exp(-1.0)

exact
Int_Trape = Trapezio(0.0,1.0,eNes,Decay)

error_T = np.abs((Int_Trape-exact)/(exact))


#%%

def PesoS(i, h, ene):
    if (i == 1 or i==ene):
        w = h/3.0
    elif(i%2 == 0):
        w = 4.0*h/3.0
    else:
        w = 2.0*h/3.0
    return w

def Simpson(a,b, eNes, fun):

    integral = []
    for ene in eNes:
        suma = 0.0
        ene_imp = ene + 1
        h = (b-a)/(ene_imp-1)
        for i in range (1,int(ene_imp)+1):
            suma += PesoS(i,h, ene_imp)*fun(a + (i-1)*h)
        integral.append(suma)
    return np.array(integral)

Int_Simp = Simpson(0,1,eNes, Decay)

error_S = np.abs((Int_Simp-exact))/(exact)

#%%

def gauss(Npoints,a,b):
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
        x[i] = x[i]*(b - a)/2.0 + (b + a)/2.0
        w[i] = w[i]*(b - a)/2.0
    return x, w
def GaussInt(a, b, eNes, fun):
    integral = []
    for ene in eNes:
        quadra = 0
        x, w = gauss(ene+1, a, b)
        for i in range(0,ene+1):
            quadra += fun(x[i])*w[i]
        integral.append(quadra)
    return np.array(integral)

IntGauss = GaussInt(0.0,1.0,eNes,Decay)

err_Gauss = np.abs(IntGauss-exact)/exact

plt.loglog(eNes, error_T, label="Trapezio")
plt.loglog(eNes, error_S, label="Simpson")
plt.loglog(eNes, err_Gauss, label="Gauss")
plt.legend()
plt.show()
