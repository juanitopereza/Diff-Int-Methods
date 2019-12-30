# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from scipy import special

#%%
a = 0.0
b = 10.0
c = 0.0
d = 10.0


# funcion a integrar
def function(x,y):
    return (4/np.pi)*np.exp(-x**2 - y**2)

# Calculamos el error absoluto entre el valor real de la integral y el calculado por los diferentes métodos.
def error(real, calculated):
    return abs((real-calculated)/real)


# Calcular la integral usando MonteCarlo
def montecarlo(f, a, b, c, d, n):
    x = np.random.random(n)*10
    y = np.random.random(n)*10
    f_eval = f(x,y)
    result = (b-a)*(d-c)*f_eval.mean()
    return np.log10(n), result

# Calcular la integral usando metodo de simpson
def trapezoide(f, a, b, c, d, n):
    x = np.linspace(a, b, n+1)
    y = np.linspace(c, d, n+1)
    h = (b-a)/n
    k = (d-c)/n
    result = 0
    for i in range(n+1):
        result += 2*(f(x[i], c) + f(x[i], d) + f(a, y[i]) + f(b, y[i]))
        for j in range(n+1):
            result += 4*(f(x[i], y[j]))
    return np.log10(n), (1/4.0)*h*k*result

# Vamos a calcular la integral
def integrate(func, a, b, c, d, n):
    # Llamamos todas la funciones que creamos para calcular la integral. Recordemos que estas funciones retornan dos numeros: el primero es el logaritmo del numero entero usado para integrar y segundo es el valor de la integral calculado.
    montecarlo_results = montecarlo(func, a, b, c, d, n)
    trapezoide_results = trapezoide(func, a, b, c, d, n)

    # creamos una lista con los resultados de las 4 funciones
    results = [montecarlo_results, trapezoide_results]
    #en la lista "numbers" guardamos los resultados de las integrales por los diferentes metodos
    numbers = [item[0] for item in results]
    #en la lista "errors" guaramos el error de nuestro cálculo con respecto al valor exacto de la integral.
    errors = [np.log10(error(exact_value, item[1])) for item in results]


    return numbers, errors
#%%

#este será el calor exacto de la integral
exact_value = (special.erf(10)-special.erf(0))**2

#creamos una lista de 20 puntos entre 0 y 10^6
points = 20
numbers = np.logspace(0, 3, points)


#estas listas vacias las usaremos para guardar nuestros valores finales
values_buffer = [[],[]]
numbers_buffer = [[],[]]

#estos nombres los usaremos para etiquetar las graficas
names = ["MonteCarlo","TRapezoide"]

#Ahora inplementaremos la función "integrate" 20 veces recorriendo 6 órdenes de magnitud.
for n in numbers:
    numbers, results = integrate(function, a, b, c, d, int(n))
    for (values_list, numbers_list, number, result) in zip(values_buffer, numbers_buffer, numbers, results): #Explicar cómo funciona "for" "in zip".
        values_list.append(result)
        numbers_list.append(number)

# Graficamos
for (values_list, numbers_list, name) in zip(values_buffer, numbers_buffer, names):
    plt.plot(numbers_list, values_list, "-o", label=name)

plt.legend(loc=3)
plt.xlabel("$\log_{10}n$")
plt.ylabel("$\log_{10}\epsilon$")
plt.grid()
plt.savefig("plots.pdf")
