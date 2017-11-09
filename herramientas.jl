
module herramientas

export Metodo_Newton 

"""Documentación Método_Newton. La función Método_Newton calcula las aproximaciones a las raíces de una función f. Sus entradas son f, la derivada de f y el punto inicial x0; se realizarán 200 iteraciones del método"""

function Método_Newton(f,df,x0)   
    x = x0                    
    for i in 1:200                   # Se crea un ciclo for para realizar las iteraciones del método de Newton
        x = x-f(x)/df(x)    
    end
    return x                         # La función arroja el valor de la aproximación de la raíz 
end

export metodoNewton2 

"""Documentación metodoNewton2. La función metodoNewton2 calcula las aproximaciones a las raíces de una función f, utilizando un ciclo while, que se ejecuta hasta que la diferencia entre los valores de la raíz calculada sea menor que una epsilon propuesta"""

function metodoNewton2(f,df,x0,epsilon) 
    xold = x0                         # Establecemos nuestra xold igual al punto inicial que nos proporcionan y posteriormente asignamos a xnew el valor de xold.
    xnew = xold
    diferencia = 10                   # Establecemos una diferencia inicial mayor que la epsilon propuesta, para que pueda entrar al ciclo e inicializamos niteracion como cero.
    niteracion = 0
    while diferencia > epsilon        # Mientras la diferencia sea mayor a epsilon, se realiza el método de Newton.
        xold = xnew                   # A x vieja se le asignará el valor de la aproximación nueva conforme se realicen las iteraciones.
        xnew = xold-f(xold)/df(xold)
        diferencia = abs(xnew-xold)   # Vemos cual es la diferencia entre la aproximación recién encontrada y la anterior.
        niteracion = niteracion + 1   # Con esta asignación, el ciclo va ir avanzando.
    end
    return xnew               # La función regresará el valor de la aproximación y el resultado de evaluar ese número en la función.
end

export metodoconintervalos 

"""Documentación metodoconintervalos. La función metodoconintervalos calcula las aproximaciones a las raíces de una función f, utilizando el método de Newton y un amplio rango de condiciones iniciales. Sus entradas son f, df y el intervalo para el rango de condiciones iniciales"""

function metodoconintervalos(f,df,intervalo)
    list=zeros(length(intervalo))  # Creamos una lista que tenga igual número de entradas que el intervalo creado con linspace.
    x = 0                          # Inicializamos
    for i in 1:length(intervalo)   # El primer ciclo va a avanzar sobre la longitud del intervalo.
        x = intervalo[i]           # A x se le va a ir asignando el valor i-ésimo del intervalo
        for n in 1:200             # Con este ciclo se va a realizar la iteración del método de Newton, con un total de 200 iteraciones.
            x = x-(f(x)/df(x))
        end
    list[i]=x;                     # En esta lista se van a ir guardando los valores de la raíz obtenida para cada iteración.
    end
    return list                           # Finalmente la función va a regresar la lista con las raíces obtenidas.
end;

export metodoNewton_epsilon

"""Documentación metodoNewton_epsilon. La función metodoNewton_epsilon calcula las aproximaciones a las raíces de una función f, utilizando el método de Newton y un amplio rango de condiciones iniciales. La función arroja la lista de las raíces que son diferentes hasta un epsilon propuesto. Sus entradas son f, df y el intervalo para el rango de condiciones iniciales."""

function metodoNewton_epsilon(f,df,intervalo)
    epsilon = 0.000000001
    t = []                                           #Arreglo vacío que almacenará las raíces que son distintas hasta un cierto epsilon.
    lista = metodoconintervalos(f,df,intervalo)      #La lista va a contener las raíces obtenidas con el método llamado metodoconintervalos, considerando el intervalo creado con linspace.
    push!(t,lista[1])                                #El primer elemento de la lista es el primer componente de la lista obtenida con el método de Newton para intervalos.
    for i in 1:length(t)                             #Con el primer ciclo for, se va a realizar el procedimiento para cada elemento de t.
        for k in 1:length(lista)                     #Con este ciclo se van a ir comparando elemento a elemento de la lista, para no incluir en el arreglo t, aquellos valores que sean iguales.
            if abs(t[i]-lista[k])>epsilon            #Si el valor absoluto de la diferencia entre las raíces es menor que el epsilon propuesta, se anexa esa raíz(componente k-ésimo de la lista) al arreglo t utilizando push
                push!(t,lista[k])
            end
        end
        return t                                     #La función nos regresará el arreglo t con las raíces que son distintas hasta un cierto epsilon.                       
    end
end

export metodo_Riemann

"""Documentación Método_Riemann. Método que permite calcular la aproximación a una integral definida sobre un intervalo [a,b] utilizando diferencias finitas. Sus entradas son la función a integrar, los extremos del intervalo y el número de particiones del mismo"""

function metodo_Riemann(f,a,b,n)   #Función que implementa el método del rectángulo, sus entradas son la función a integrar, los extremos del intervalo y el número de segmentos en que se va a dividir el intervalo original.
    intervalo=(b-a)/n              #Considerando que los segmentos tengan la misma longitud, podemos obtener esta longitud haciendo el cociente de la resta de los extremos del intervalo entre el número de segmentos en que queremos separar el intervalo.
    suma=0                         #En suma se va a ir guardando el valor de la implementación del método del rectángulo para cada segmento del intervalo. Inicializamos suma como cero.
    x=zeros(1,n)                   #Aquí se van a ir guardando los valores de los puntos de cada segmento.
    x[1]=a                         #El primer punto del intervalo es el extremo a
    for i in 2:n                   #Con este ciclo for vamos a ir avanzando sobre los puntos del intervalo
        x[i]=x[i-1]+intervalo      #El punto nuevo del intervalo se va a obtener sumando el punto anterior del intervalo (que para la primera vuelta del ciclo es a) más la longitud del segmento.
        suma = suma + (x[i]-x[i-1])*f((x[i-1]+x[i])/2)   #Aquí vamos a anexar a suma el valor obtenindo aplicando la fórmula mostrada arriba para cada segmento del intervalo.
    end
    return suma                    #Finalmente la función nos va a arrojar el valor de la suma, que corresponde a la aproximación de la integral utilizando el método del rectángulo.
end

export metodo_del_trapecio

"""Documentación Método_del_trapecio. Método que permite calcular la aproximación a una integral definida sobre un intervalo [a,b] utilizando el método del trapecio. Sus entradas son la función a integrar, los extremos del intervalo y el número de particiones del mismo"""

function metodo_del_trapecio(f,a,b,n)   #Función que implementa el método del trapecio, sus entradas son la función a integrar, los extremos del intervalo y el número de segmentos en que se va a dividir el intervalo original.
    intervalo=(b-a)/n              #Considerando que los segmentos tengan la misma longitud, podemos obtener esta longitud haciendo el cociente de la resta de los extremos del intervalo entre el número de segmentos en que queremos separar el intervalo.
    suma=0                         #En suma se va a ir guardando el valor de la implementación del método del del trapecio para cada segmento del intervalo. Inicializamos suma como cero.
    x=zeros(1,n)                   #Aquí se van a ir guardando los valores de cada segmento.
    x[1]=a                         #El primer punto del intervalo es el extremo a
    for i in 2:n                   #Con este ciclo for vamos a ir avanzando sobre los puntos del intervalo
        x[i]=x[i-1]+intervalo      #El punto nuevo del intervalo se va a obtener sumando el punto anterior del intervalo (que para la primera vuelta del ciclo es a) más la longitud del segmento.
        suma = suma + ((x[i]-x[i-1]))*((f(x[i-1])+f(x[i]))/2)   #Aquí vamos a anexar a suma el valor obtenindo aplicando la fórmula mostrada arriba para cada segmento del intervalo.
    end
    return suma                    #Finalmente la función nos va a arrojar el valor de la suma, que corresponde a la aproximación de la integral utilizando el método del trapecio.
end

export metodo_de_Simpson

"""Documentación Método_de_Simpson. Método que permite calcular la aproximación a una integral definida sobre un intervalo [a,b] utilizando el método de Simpson. Sus entradas son la función a integrar, los extremos del intervalo y el número de particiones del mismo"""

function metodo_de_Simpson(f,a,b,n) #Función que implementa el método de Simpson, sus entradas son la función a integrar, los extremos del intervalo y el número de segmentos en que se va a dividir el intervalo original.
    intervalo=(b-a)/n               #Considerando que los segmentos tengan la misma longitud, podemos obtener esta longitud haciendo el cociente de la resta de los extremos del intervalo entre el número de segmentos en que queremos separar el intervalo.
    suma=0                           #En suma se va a ir guardando el valor de la implementación del método del Simpson para cada segmento del intervalo. Inicializamos suma como cero.
    x=zeros(1,n)                     #Aquí se van a ir guardando los valores de los puntos para cada segmento.
    x[1]=a                           #El primer punto del intervalo es el extremo a
    for i in 2:n                     #Con este ciclo for vamos a ir avanzando sobre los puntos del intervalo
        x[i]=x[i-1]+intervalo        #El punto nuevo del intervalo se va a obtener sumando el punto anterior del intervalo (que para la primera vuelta del ciclo es a) más la longitud del segmento.
        suma = suma + (1/6)*(x[i]-x[i-1])*(f(x[i-1])+4*f((x[i-1]+x[i])/2)+f(x[i]))    #Aquí vamos a anexar a suma el valor obtenindo aplicando la fórmula mostrada arriba para cada segmento del intervalo.
    end
    return suma                       #Finalmente la función nos va a arrojar el valor de la suma, que corresponde a la aproximación de la integral utilizando el método de Simpson.
end

export Método_Euler_Explícito

"""Documentación del método de Euler explícito. Método de integración numérica que permite resolver ecuaciones diferenciales ordinarias. Sus entradas son la función f, los tiempos inicial t0 y final tf, el tamaño de paso h y el punto inicial x0."""

function Método_Euler_Explícito(f,t0,x0,tf,h)  #Función eulerMethod cuyas entradas son la función f, los tiempos inicial y final, el punto inicial x0 y el tamaño de paso h.
    t = linspace(t0,tf,(tf-t0)/h)   #Creamos un intervalo cuyos extremos son los tiempos inicial y final, que tiene n elementos, considerando n=(tf-t0)/h
    x=zeros((tf-t0)/h)            #Creamos un arreglo de n elementos donde se irán guardando las x´s obtenidas con el método de Euler.
    x[1]=x0                         #Condición inicial
    for i in 2:length(t)            #Ciclo for para implementar la fórmula de recurrencia mostrada arriba.
        x[i] = x[i-1] + h*f(t[i-1],x[i-1])
    end
    return(t,x)                     #La función eulerMethod arroja los valores obtenidospara t y x.
end

export Metodo_Euler_explícito_dimensiones

"""Documentación del método de Euler explícito, independiente de las dimensiones del sistema. Sus entradas son la función f, el punto inicial x0 y list que contiene los tiempos inicial y final y el tamaño de paso h."""

function Metodo_Euler_explícito_dimensiones(f,list,x0)    # Función metodo_euler cuyas entradas son la función f, list que es el intervalo de tiempo y el punto inicial x0.
    x = x0                          # Condición inicial
    h = list[2]-list[1]             # El tamaño de paso h se obtiene de restar los primeros 2 elementos de la lista.
    listx = []                      # Arreglo vacío que guardará los resultados obtenidos con este método.
    push!(listx,x)                  # Agregamos al arreglo listx, el punto inicial con push
    for i in 2:length(list)         # Ciclo for que implementará la fórmula de recurrencia mostrada arriba para el método de Euler 
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x)              # Para cada vuelta del ciclo, se agrega el resultado obtenido al arreglo listx 
    end
    return listx                    # La función metodo_euler regresa listx
end

export Método_Euler_modificado

"""Documentación del método de Euler modificado, o punto medio. Método de integración numérica que permite resolver ecuaciones diferenciales ordinarias. Sus entradas son la función f, los tiempos inicial y final, el tamaño de paso y el punto inicial x0."""

function Método_Euler_modificado(f,t0,x0,tf,h)  # Función Método_Euler_modificado cuyas entradas son la función f, los tiempos inicial y final, el punto inicial x0 y el tamaño de paso h.
    t = linspace(t0,tf,(tf-t0)/h)     # Creamos un intervalo cuyos extremos son los tiempos inicial y final, con n elementos, considerando n=(tf-t0)/h
    t1=zeros((tf-t0)/h)               # Creamos un arreglo de n elementos donde se irán guardando los tiempos medios obtenidos.
    x=zeros((tf-t0)/h)                # Creamos un arreglo de n elementos donde se irán guardando las x´s obtenidas con el método de Euler modificado.
    x[1]=x0                           # Condición inicial
    t1[1]=t0                          # El primer elemento del arreglo para los tiempos medios es igual al tiempo inicial.
    for i in 1:length(t)-1            # Ciclo for para implementar la fórmula de recurrencia mostrada arriba.
        t1[i+1]=(t[i+1]+t[i])/2       # Fórmula para calcular los tiempos medios.
        x[i+1] = x[i] + h*f(x[i]+(h/2)*f(x[i],t[i]),t1[i+1])                       
    end
    return(t1,x)                      # La función Método_Euler_modificado arroja los valores obtenidos para t1 y x.
end

export metodoNewton 

function metodoNewton(f,df,x0,t)   
    x = x0                    
    for i in 1:200                   # Se crea un ciclo for para realizar las iteraciones del método de Newton
        x = x-f(x)/df(x)    
    end
    return x                         # La función arroja el valor de la aproximación de la raíz 
end

export Método_implícito_Euler_dimensiones

"""Documentación del método de Euler implícito. Método de integración numérica que permite resolver ecuaciones diferenciales ordinarias. Sus entradas son la función f, la derivada de la función f, los tiempos inicial y final, el tamaño de paso y el punto inicial x0. Este método hace uso del método de Newton para aproximación de raíces."""

function Método_implícito_Euler(f,df,t0,tf,h,x0)  # Función Método_implícito_Euler cuyas entradas son la función f, la derivada de la función f, los tiempos inicial y final, el punto inicial x0 y el tamaño de paso h.
    t = linspace(t0,tf,(tf-t0)/h)                 # Creamos un intervalo cuyos extremos son los tiempos inicial y final, que tiene n elementos, considerando n=(tf-t0)/h
    listax=zeros((tf-t0)/h)                       # Creamos un arreglo de n elementos donde se irán guardando las x´s obtenidas con el método.
    x = x0                                        # Condición inicial
    listax[1] = x0
    for i in 2:length(t)                          # Ciclo for para implementar la fórmula de recurrencia mostrada arriba.
        g(x) = x - listax[i-1] - h*f(x,t[i])
        dg(x) = 1 - h*df(x,t[i])
        x = metodoNewton(g,dg,listax[i-1],t[i])   # Utilizando el método de Newton mostrado arriba
        listax[i] = x
    end
    return(t,listax)                              # La función Método_implícito_Euler arroja los valores obtenidos para t y listax.
end

export runge_kutta_4

"""Documentación del método de Runge-Kutta de orden 4. Método de integración numérica que permite resolver ecuaciones diferenciales ordinarias. Sus entradas son la función f, los tiempos inicial t0 y final tf, el tamaño de paso h y el punto inicial x0."""

function runge_kutta_4(f,x0,t0,tf,h)  # Función Runge-Kutta cuyas entradas son la función f, los tiempos inicial y final, el punto inicial x0 y el tamaño de paso h.
    n=round((tf-t0)/h)+1               
    listt=linspace(t0,tf,n) 
    listx = zeros(n)                  # Creamos un arreglo de n elementos, considerando n=(tf-t0)/h, donde se irán guardando las x´s obtenidas con el método de Runge Kutta.
    listx[1] = x0 
    for i in 1:length(listx)-1        # Ciclo for para implementar la fórmula de recurrencia mostrada arriba.
        k1 = f(listx[i], listt[i])
        k2 = f(listx[i] + h*(k1)/2, listt[i+1])
        k3 = f(listx[i] + h*(k2)/2, listt[i+1])
        k4 = f(listx[i] + h*(k3), listt[i],)
        listx[i+1] = listx[i] + h/6*(k1 + 2*(k2) + 2*(k3) + k4)
    end
    return listt,listx    # La función Runge-Kutta arroja los valores obtenidos para listt y listx.
end

export Metodo_Runge_Kutta_dimensiones

"""Documentación del método de Runge-Kutta de orden 4, independiente de las dimensiones del sistema. Sus entradas son la función f, el punto inicial x0 y list que contiene los tiempos inicial t0 y final tf y el tamaño de paso h."""

function Metodo_Runge_Kutta_dimensiones(f,list,x0)    # Función metodo_Runge_Kutta cuyas entradas son la función f, list que es el intervalo de tiempo y el punto inicial x0.
    x = x0                          # Condición inicial
    h = list[2]-list[1]             # El tamaño de paso h se obtiene de restar los primeros 2 elementos de la lista.
    listx = []                      # Arreglo vacío que guardará los resultados obtenidos con este método.
    push!(listx,x)                  # Agregamos al arreglo listx, el punto inicial con push
    for i in 2:length(list)         # Ciclo for que implementará las fórmulas de recurrencia mostrada arriba para el método de Runge-Kutta 
        t = i*h
        k1 = f(x,t);
        k2 = f(x+(h/2)*k1,t+(h/2));
        k3 = f(x+(h/2)*k2, t+(h/2));
        k4 = f(x+h*k3, t);
        x = x+(h/6)*(k1+2*k2+2*k3+k4);
        push!(listx,x)              # Para cada vuelta del ciclo, se agrega el resultado obtenido al arreglo listx 
    end
    return listx                    # La función metodo_euler regresa listx
end
end
