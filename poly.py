# -*- coding: utf-8 -*-
from pylab import *
from numpy import *
from sys import exit


class Poligono():

    """
    Clase que genera un poligono regular con centroide c,
    numero de lados n, y escala de lados l (radio del circuncirculo)

    El centroide hay que pasarlo como lista o array c = [a,b] / c = array([a,b])
    """

    def __init__(self, c, n, l, regular = True, vertex_list = None):

        """
        Constructor de la clase, pasa los argumentos de centroide,
        numero de lados y radio del circulo que circunscribe al poligono.

        Tambien define los vectores tangentes y normales a los lados del poligono
        """


        ######### Define Variables internas del poligono


        self.n = n

        self.vertices = []
        self.tangentes = []
        self.normales = []

        self.regular = regular  #es necesario??

        ################################################ Vertices

        if regular:

            """
            Los siguietes angulos solo se definen si se tiene un poligono regular
            """

            self.l = l
            self.centroide = array(c)
            self.angulo = ( (self.n - 2.)*pi )/ self.n #angulo interior
            self.angulo_sector = ( 2.0*pi )/ n


            if n%2 != 0:
                theta = 0.0
            else:
                theta = self.angulo_sector / 2

            for i in range(self.n):

                vertice = l * array( [cos(theta), sin(theta)] )
                self.vertices.append(vertice)
                theta += self.angulo_sector

            self.vertices = array(self.vertices)

            self.longitud_de_lado = norm( self.vertices[0] - self.vertices[1] )
            self.longitud_de_arco = n* self.longitud_de_lado

        else:
            """si el poligono no es regular: """

            for i in range(len(vertex_list)):
                self.vertices.append( array(vertex_list[i]).astype("float") )

            x = 0.0
            y = 0.0
            for vertice in self.vertices:

                x += vertice[0]
                y += vertice[1]

            x /= self.n
            y /= self.n

            self.centroide = array([x, y])


            self.l = max(map( norm, self.centroide - array(self.vertices)))


            self.longitud_de_lado = [] # Para un poligono no-regular

            for i in range(self.n):
                if i != (self.n - 1):
                    self.longitud_de_lado.append( norm(self.vertices[i] - self.vertices[i+1]) )
                else:
                    self.longitud_de_lado.append( norm(self.vertices[i] - self.vertices[0]) )

            self.longitud_de_arco = sum(self.longitud_de_lado)


            #self.angulo_sector =
            """esto se puede definir para poligonos iregulares, pero se necita una lista que corresponda a cada vertice"""

            if len(self.vertices) != self.n:
                print "\nEl numero de vertices no coincide con los lados del poligono"
                print "verificar lista de vertices e intentar de nuevo\n"
                exit(0)



        #################################################### Vectores Tangentes y normales

        for i in range(self.n):

            if i != (self.n -1):
                tang = array(array(self.vertices[i+1]) - array(self.vertices[i])).astype("float")
            else:
                tang = array(array(self.vertices[0]) - array(self.vertices[i])).astype("float")

            tang /= norm(tang)

            a, b = tang
            nor = -b, a

            self.tangentes.append(array(tang))
            self.normales.append(array(nor))

    def pos_de_arco_regular(self, x, frontera):

        """
        Regresa la posici髇 de arco del punto x de colision, conociendo la frontera con la
        cual ha colisionado si se trata de un poligono regular
        """


        s = norm( x - self.vertices[ frontera ] ) + self.longitud_de_lado*( frontera )

        return s


    def pos_de_arco_no_regular(self, x, frontera):

        """
        Regresa la posici贸n de arco del punto x de colision, conociendo la frontera con la
        cual ha colisionado si se trata de un poligono NO regular
        """
        s = norm( x - self.vertices[ frontera ] ) + sum(self.longitud_de_lado[:frontera])

        return s


    def vector_normal(self, x, frontera):

        """
        Devuelve el vector normal al poligon en el punto x
        """

        v_normal = self.normales[ frontera ]

        return v_normal

    def vector_tangente(self, x, frontera):

        """
        Devuelve el vector tangente al poligono en el punto x
        """

        v_tangente = self.tangentes[ frontera ]

        return v_tangente

    def draw(self):

        """
        Dibuja en matplotlib el poligono
        """

        axis("equal")

        ls = list(self.vertices)
        ls.append(self.vertices[0])
        ls = array(ls)

        plot(ls[:,0], ls[:,1], "-")

        #show()


    def tiempo_interseccion(self, r, v, Orbita):


        """
        Calcula el tiempo de interseccion de una particula con pos r y velocidad v
        con el poligono.

        Devuelve, el timepo de interseccion y la frontera con la qua acaba de ocurrir una colision!!!
        """

        tiempos_posibles = []
        lista_fronteras = []



        for i in range(self.n):

            if i == Orbita.frontera_previa:
                continue  # evita colisiones consecutivas con la misma frontera

            p = self.vertices[i]  # punto sobre recta (el que sea, pero uso uno de los vertices por facilidad)
            n = self.normales[i]  # aqui la n se refiere a un vector normal (NO al numero de lados !!!)

            t = dot( (p - r), n) / dot( v, n ) #formula para tiempo de interseccion

            if sign(t) >= 0:

                tiempos_posibles.append(t)
                lista_fronteras.append(i)

        tiempo_int = min(tiempos_posibles)
        numero_en_lista = argmin(tiempos_posibles)
        frontera_previa = lista_fronteras[numero_en_lista]

        return tiempo_int, frontera_previa



    def pos_ini(self):

        """
        el poligono para iniciar la dinamica
        """

        #x0 = ((rand(2) + self.centroide ) - array([0.5, 0.5]))*2*self.l
        x0 = ( rand(2) - array([0.5, 0.5]) )*2*self.l + self.centroide #+ self.centroide

        #x0 = rand(2)*2*self.l + self.centroide

        while True:

            L1 = [] # lista de vectores que van de los vertices al punto x0
            L2 = [] # proyecciones de vect de vert a x0 sobre normal correspondiente
            j = 0

            for i in range(self.n):
                a = x0 - self.vertices[i]
                L1.append(a)

                b = dot(a, self.normales[i])
                L2.append(b)

            L2 = array(L2)

            if all(L2 > 0):
                break

            x0 = ( rand(2) - array([0.5, 0.5]) )*2*self.l + self.centroide
            #plot(x0[0], x0[1], "o")

        return x0


    def vel_ini(self):

        """
        Aleatoriamente da una velocidad inicial unitaria
        """
        ang = rand()*2*pi
        v0 = array([cos(ang), sin(ang)])
        return v0


class ReglasDeColision():

    """
    Clase que incorpora las reglas de colision posibles.
    Toman los vectores normal, tangente y la velocidad y devuelven una velocidad nueva NORMALIZADA.

    La sintaxis de todas las funciones en esta clase deberia de ser la misma para reducir posibles errores
    """


    def david(self, x, v, frontera, mesa, alfa):

        """
        Regla de colision de David, re-escala el angulo de reflexion por el factor alfa
        x - pt de colision
        v - vel de colision
        frontera - se refiere al lado con el que se acaba de colisionar...
        mesa - mesa donde se lleva a cabo la dinamica
        alfa - coef de inelasticidad
        """

        nor = mesa.vector_normal(x, frontera)
        tan = mesa.vector_tangente(x, frontera)

        proy_tangente = dot(tan, v)
        proy_normal = dot(nor, v)

        signo = sign(proy_tangente)

        angulo_entrada = arccos( dot(-v, nor) )# Es -v, porque la convencion que estoy usando es que las normales apuntan hacia adentro
        self.angulo_entrada = angulo_entrada

        angulo_salida = signo*angulo_entrada * alfa
        self.angulo_salida = angulo_salida

        componente_normal = cos(angulo_salida)
        componente_tangente = sin(angulo_salida)

        v_nuevo = (nor*componente_normal) + (tan*componente_tangente)
        v_nuevo /= norm(v_nuevo)

        return v_nuevo


    def alonso(self, x, v, frontera, mesa, alfa):

        """
        Mi regla de colision, re-escala la componente tangente de la velocidad
        por el factor alfa y luego la normaliza
        """

        nor = mesa.vector_normal(x, frontera)
        tan = mesa.vector_tangente(x, frontera)

        proy_tangente = dot(tan, v)
        proy_normal = dot(nor, v)

        signo = sign(proy_tangente)

        angulo_entrada = arccos( dot(-v, nor) )# Es -v, porque la convencion que estoy usando es que las normales apuntan hacia adentro
        self.angulo_entrada = angulo_entrada

        v_nuevo = ( proy_tangente * alfa * tan) - ( proy_normal* nor)

        v_nuevo /= norm(v_nuevo)

        angulo_salida = arccos( dot(v_nuevo, nor) )  # esta haciendo algo muy raro...

        self.angulo_salida = signo*angulo_salida

        #print dot(v_nuevo, nor)
        #print self.angulo_salida

        return v_nuevo



    def __init__(self, regla):

        """Inicializa el objeto "extrano" ReglasDeColision, que contiene las reglas de colision que se van a utilizar
        y debido a que en la funcion de reflexion se tiene que calcular los angulos de entrada y salida que seran necesarios
        despues, se inicializan estas variables como 0
        """

        self.angulo_entrada = 0.0
        self.angulo_salida = 0.0

        if regla == "d":
            self.colision = self.david
        elif regla == "a":
            self.colision = self.alonso





class MapeosTangentes():


    #def A_a(self,n, p, alfa, theta):
        #"""
        #Devuelve la derivada del mapeo de colisi贸n para mi regla
        #"""

        #np = n*p
        #matriz_aux = array([n[0], n[1], p[1], p[0]]).reshape(2,2)

        #dbeta_dp = []
        #for i in range(2):
            #j = 1 - i
            #dbeta_dp.append(n[i]*((alfa**2.0)*(n[i]*p[i]-n[j]*p[j]) - dot(n,p))/ ( (alfa**2)*det(matriz_aux) + dot(n,p)**2)**(3.0/2.0) )

        #M_A = array([1.0,0.0,0.0,0.0,  0.0, 1.0, 0.0, 0.0,   0.0,0.0,(alfa*(n[1])**2 - n[0]**2)*dbeta_dp[0], -dot(n,n)*dbeta_dp[1],   0.0,0.0,-dot(n,n)*dbeta_dp[0], (alfa*(n[0])**2 - n[1]**2)*dbeta_dp[1]]).reshape(4,4)

        #return M_A


    #def A_d(self,n, p, alfa, theta):
        #"""
        #Devuelve la derivada del mapeo de colisi贸n para la regla de David, theta deber谩 ser el 谩ngulo de salida ya modificado
        #"""

        #raiz = sqrt(1 - dot(p,n)**2)

        #a = []
        #b = []

        #for i in range(2):

            #a.append( alfa*n[i]*(n[1]*cos(theta) - n[0]*sin(theta))/raiz )
            #b.append( -alfa*n[i]*(n[0]*cos(theta) + n[1]*sin(theta) )/raiz )

        #M_D = array([1.0,0.0,0.0,0.0,  0.0, 1.0, 0.0, 0.0,    0.0,0.0,a[0], a[1],     0.0, 0.0, b[0], b[1]]).reshape(4,4)

        #return M_D

    #def m_l(self,t):
        #"""
        #Devuelve la matriz correspondiente a la transf. de los vectores tangentes
        #durante el tiempo de vuelo libre t
        #"""

        #matriz = array([1.0, 0.0, t, 0.0,    0.0, 1.0, 0.0, t,    0.0,0.0,1.0,0.0,   0.0,0.0,0.0,1.0]).reshape(4,4)

        #return matriz


    def Df1(self, n, v):

        n1, n2 = n

        matriz = array([1-2*n1**2, -2*n1*n2, -2*n1*n2, 1-2*n2**2]).reshape(2,2)

        return matriz


    def Df2(self, v):

        matriz = array([1,0,0,1,1,0,0,1]).reshape(4,2)

        return matriz

    def Df3(self, v):

        v1, v2, v3, v4 = v

        d_theta_1 = -v2/sqrt(v1**2 + v2**2)
        d_theta_2 = v1/sqrt(v1**2 + v2**2)

        matriz = array([ -d_theta_1, -d_theta_2, 0, 0,0,0,1,0,0,0,0,1]).reshape(3,4)

        return matriz

    def Df4(self, alfa, v):

        matriz = array([(1-alfa), 0, 0, 0, 1, 0, 0, 0, 1]).reshape(3,3)

        return matriz

    def Df5(self, v):

        v1, v2, v3 = v

        matriz = array([ -sin(v1), 0,0, -cos(v1), 0,0, cos(v1),0,0, -sin(v1),0,0,0,1,0,0,0,1]).reshape(6,3)

        return matriz

    def Df6(self, v):

        v1, v2, v3, v4, v5, v6 = v

        matriz = array([ v5, v6, 0, 0, v1, v2, 0, 0, v5, v6, v3, v4]).reshape(2, 6)

        return matriz

    def transformacion_colision_david(self, mesa, frontera, time, dG, normas, p_antes, p_despues, theta, alfa):
        """
        Trasforma los vectores tangentes de acuerdo a la regla de colision de David y agrega a la lista que se le pase como
        normas la norma después del proceso Gram-Shmidt. t = tiempo de vuelo libre
        """

        t = mesa.tangentes[frontera]
        n = mesa.normales[frontera]


        #dot_pn = dot(p_antes,n)
        #dot_pt = dot(p_antes, tangent_vector)


        p = p_antes

        y1 = p - 2.0*dot(p,n)*n

        y2 = array([y1[0], y1[1], y1[0], y1[1]])

        y3 = array([arctan(dot(y2[:2],t)/dot(y2[:2],n)), y2[0], y2[1]])

        y4 = array([(1- alfa)*y3[0], y3[1], y3[2] ])

        y5 = array([ cos(y4[0]), -sin(y4[0]), sin(y4[0]), cos(y4[0]), y4[1], y4[2]  ])

        y6 = dot(y5[:4].reshape(2,2) , y5[4:])

        #print "y6: ", y6, "p_despues: ", p_despues
        #print "diferencia: ", y6 - p_despues
        #print  "angulo: ", theta, "y4[0]: ", y4[0]


        #print "\n",theta
        #theta = -arctan(dot_pt/dot_pn)
        #print theta
        #print
        #c = cos(-theta*(1-alfa))
        #s = sin(-theta*(1-alfa))



        ##derivadas de las funciones sin(theta), cos(theta) con respecto a p
        #diff_c_p1 = -s*(1 - alfa)*(-p2/norm(p)**2)
        #diff_c_p2 = -s*(1 - alfa)*(p1/norm(p)**2)

        #diff_s_p1 = c*(1 - alfa)*(p2/norm(p)**2)
        #diff_s_p2 = c*(1 - alfa)*(p1/norm(p)**2)


        ##derivadas que conforman la derivada total de la componente inelastica del mapeo
        #diff_f1_p1 = c + p1*diff_c_p1 - p2*diff_s_p1
        #diff_f1_p2 = -s + p1*diff_c_p2 - p2*diff_s_p2

        #diff_f2_p1 = s + p1*diff_s_p1 + p2*diff_c_p1
        #diff_f2_p2 = c + p1*diff_s_p2 + p2*diff_c_p2



        #F_inel = array([diff_f1_p1, diff_f1_p2, diff_f2_p1, diff_f2_p2]).reshape(2,2)

        #print dot_pn, dot_pt, theta, c, s


        #df1dp1 = (n[1]**2 -n[0]**2)*c + 2*n[0]*n[1]*s - p[0]*(n[1]**2 -n[0]**2)*(1-alfa)*s*p[1] - p[0]*2*n[0]*n[1]*(1-alfa)*c*p[1] - p[1]*(n[1]**2 -n[0]**2)*c*(1-alfa)*p[1] + p[1]*2*n[0]*n[1]*s*(1-alfa)*p[1]# - ((n[1]**2 -n[0]**2)*s -2*n[0]*n[1]*c)*(p[0]/sqrt(1-p[0]**2))

        #df1dp2 = (n[1]**2 -n[0]**2)*s - 2*n[0]*n[1]*c + p[1]*(n[1]**2 -n[0]**2)*(1-alfa)*c*p[0] - p[1]*2*n[0]*n[1]*(1-alfa)*s*p[0] + p[0]*(n[1]**2 -n[0]**2)*s*(1-alfa)*p[0] + p[0]*2*n[0]*n[1]*c*(1-alfa)*p[0]# - ((n[1]**2 -n[0]**2)*c +2*n[0]*n[1]*s)*(p[1]/sqrt(1-p[1]**2))

        #df2dp1 = (n[1]**2 -n[0]**2)*s - 2*n[0]*n[1]*c - p[0]*(n[1]**2 -n[0]**2)*(1-alfa)*c*p[1] + p[0]*2*n[0]*n[1]*(1-alfa)*s*p[1] - p[1]*(n[0]**2 -n[1]**2)*s*(1-alfa)*p[1] + p[1]*2*n[0]*n[1]*c*(1-alfa)*p[1]# - ((n[0]**2 -n[1]**2)*c - 2*n[0]*n[1]*s)*(p[0]/sqrt(1-p[0]**2))

        #df2dp2 = (n[0]**2 -n[1]**2)*c - 2*n[0]*n[1]*s + p[1]*(n[0]**2 -n[1]**2)*(1-alfa)*s*p[0] - p[1]*2*n[0]*n[1]*(1-alfa)*c*p[0] + p[0]*(n[1]**2 -n[0]**2)*c*(1-alfa)*p[0] - p[0]*2*n[0]*n[1]*s*(1-alfa)*p[0]# - ((n[1]**2 -n[0]**2)*s - 2*n[0]*n[1]*c)*(p[1]/sqrt(1-p[1]**2))

        #df2dp2 = (n[0]**2 -n[1]**2)*c + 2*n[0]*n[1]*s + p[1]*(n[0]**2 -n[1]**2)*(n[1]+n[0])*(1-alfa)*s + 2*p[1]*n[0]*n[1]*(n[1]+n[0])*(1-alfa)*c + p[0]*(n[1]**2 -n[0]**2)*c*(1-alfa)*(n[1]+n[0]) + 2*p[0]*n[0]*n[1]*s*(1-alfa)*(n[1]+n[0]) - ((n[1]**2 -n[0]**2)*s - 2*n[0]*n[1]*c)*(p[1]/sqrt(1-p[1]**2))



        #A_d = array([df1dp1, df1dp2, df2dp1, df2dp2]).reshape(2,2)
        #dot_pt = dot(p_antes,t)

        #theta = -arctan(dot_pt / dot_pn)


        #d_theta = n*dot_pt - t*dot_pn      #cada termino es el factor de dp_i pues d_theta_real = sum_i (n_i(p.t) - t_i(p.n))

        #da = alfa*(n[0]*cos(alfa*theta) - n[1]*cos(alfa*theta)) * d_theta
        #db = alfa*(n[0]*cos(alfa*theta) + n[1]*sin(alfa*theta)) * d_theta

        #print da
        #print db


        #raiz = sqrt(1 - dot(p_antes,n)**2)

        #sgn = sign(dot(p_antes, t))
        #print "El signo es:\t", sgn

        #a = []
        #b = []

        #for i in range(2):

            #a.append( (alfa*n[i]*(n[1]*sgn*cos(theta) - n[0]*sin(theta)))/ raiz )
            #b.append(-(alfa*n[i]*(n[0]*sgn*cos(theta) + n[1]*sin(theta)))/ raiz )

        #A_d = array([da[0], da[1], db[0], db[1]]).reshape(2,2)

        #print "a:\t", a[0], "\t", a[1]
        #print "b:\t", b[0], "\t", b[1]

        delta_p = p_despues - p_antes


        dG_nuevo = []

        for dg in dG:

            dq1, dq2, dp1, dp2 = dg
            dq = array([dq1, dq2])
            dp = array([dp1, dp2])

            dot_nq = dot(n, dq)

            dq1_nuevo , dq2_nuevo = dq +  delta_p*dot_nq/dot(p_antes, n) # 2.0*dot_nq*n

            #dq1_nuevo = dq1 - 2.0*dot_nq*n[0]
            #dq2_nuevo = dq2 - 2.0*dot_nq*n[1]

            #dp1_nuevo = a[0]*dp1 + a[1]*dp2
            #dp2_nuevo = b[0]*dp1 + b[1]*dp2

            #A_d = array([a[0], a[1], b[0], b[1]]).reshape(2,2)

            #dp_nuevo = dot(A_d, dp)
            #dp1_nuevo, dp2_nuevo = dp_nuevo
            #dp1_nuevo = df1dp1 * dp1 + df1dp2 * dp2
            #dp2_nuevo = df2dp1 * dp1 + df2dp2 * dp2

            #print "D dp1:\t", dp1 - dp1_nuevo
            #print "D dp2:\t", dp2 -dp2_nuevo



            #dp_nuevo = dp - 2.0*dot(dp, n)*n #Parte elástica del mapeo
            #dp1_nuevo, dp2_nuevo = dp_nuevo

            ##derivadas de las funciones sin(theta), cos(theta) con respecto a p
            #diff_c_p1 = -(-s*(1.0 - alfa))*(-dp2_nuevo/norm(dp_nuevo)**2)
            #diff_c_p2 = -(-s*(1.0 - alfa))*(dp1_nuevo/norm(dp_nuevo)**2)

            #diff_s_p1 = -(c*(1.0 - alfa))*(-dp2_nuevo/norm(dp_nuevo)**2)
            #diff_s_p2 = -(c*(1.0 - alfa))*(dp1_nuevo/norm(dp_nuevo)**2)


            ##derivadas que conforman la derivada total de la componente inelastica del mapeo ya evaluadas en f_el(dp)...
            #diff_f1_p1 = c + dp1_nuevo*diff_c_p1 - dp2_nuevo*diff_s_p1
            #diff_f1_p2 = -s + dp1_nuevo*diff_c_p2 - dp2_nuevo*diff_s_p2

            #diff_f2_p1 = s + dp1_nuevo*diff_s_p1 + dp2_nuevo*diff_c_p1
            #diff_f2_p2 = c + dp1_nuevo*diff_s_p2 + dp2_nuevo*diff_c_p2



            #F_inel = array([diff_f1_p1, diff_f1_p2, diff_f2_p1, diff_f2_p2]).reshape(2,2)


            #dp_nuevo = dot(F_inel, dp_nuevo) #parte inelastica del mapeo

            #dp1_nuevo, dp2_nuevo = dp_nuevo

            #dg_nuevo = array([dq1_nuevo, dq2_nuevo, dp1_nuevo, dp2_nuevo])

            #dp_nuevo = dot( self.Df6(y5), dot(self.Df5(y4), dot( self.Df4(y3), dot( self.Df3(y2), dot( self.Df2(y1), dot( self.Df1(n,p), dp ))))))

            m = dot(self.Df6(y5), self.Df5(y4))
            m = dot(m, self.Df4(alfa, y3))
            m = dot(m, self.Df3(y2))
            m = dot(m, self.Df2(y1))
            m = dot(m, self.Df1(n,p))

            dp_nuevo = dot(m, dp)

            dp1_nuevo, dp2_nuevo = dp_nuevo

            dg_nuevo = array([dq1_nuevo, dq2_nuevo, dp1_nuevo, dp2_nuevo])

            dG_nuevo.append(dg_nuevo)

        dG = dG_nuevo

        return array(dG)

    #def transformacion_colision_alonso(self, mesa, frontera, t, dG, normas, p_actual, p_anterior, theta, alfa):
        #"""
        #Trasforma los vectores tangentes de acuerdo a la regla de colision de David y agrega a la lista que se le pase como
        #normas la norma despu茅s del proceso Gram-Shmidt
        #"""

        #n = mesa.normales[frontera]

        #np = n*p

        #dot_nn = dot(n,n)

        #dot_np = dot(n,p)

        #matriz_aux = array([n[0], n[1], p[1], p[0]]).reshape(2,2)

        #dbeta_dp = []
        #for i in range(2):
            #j = 1 - i
            #dbeta_dp.append(n[i]*((alfa**2.0)*(n[i]*p[i]-n[j]*p[j]) - dot_np)/ ( (alfa**2)*det(matriz_aux) + dot_np**2)**(3.0/2.0) )
            ##estas son las derivadas de beta con respecto a p_i

        #for dg in dG:

            #dq1, dq2, dp1, dp2 = dg
            #dq = array([dq1, dq2])

            #dot_nq = dot(n, dq)

            #dq1_nuevo , dq2_nuevo = dq - 2.0*dot_nq*n

            #dp1_nuevo = ((alfa*(n[1])**2 - n[0]**2)*dbeta_dp[0])*dp1 + (-dot_nn*dbeta_dp[1])*dp2

            #dp2_nuevo = (-dot_nn*dbeta_dp[0])*dp1 + ((alfa*(n[0])**2 - n[1]**2)*dbeta_dp[1])*dp2

            #dg = array([dq1_nuevo, dq2_nuevo, dp1_nuevo, dp2_nuevo])

        #return dG




    def __init__(self, regla):

        """
        Elige la función colision corrrespondiente a la regla elegida.

        Se define aqui el constructor para poder elegir entre las funciones ya definidas
        anteriormente.
        """

        if regla == "d":
            self.transformacion_colision = self.transformacion_colision_david
        elif regla == "a":
            self.transformacion_colision = self.transformacion_colision_alonso




    def evolucion_tangentes(self, mesa, frontera, t, dG, normas, v_antes, v_despues, theta, alfa):

        """
        Evoluciona los vectores tangentes en el tiempo de vuelo libre y
        debidoa a la colision
        """

        #print "dG antes:   ", list(dG)
        dG_nuevo = []
        for dg in dG:
            dq1, dq2, dp1, dp2 = dg

            dq1_nuevo = dq1 + t*dp1
            dq2_nuevo = dq2 + t*dp2

            dg_nuevo = [dq1_nuevo, dq2_nuevo, dp1, dp2]

            dG_nuevo.append(dg_nuevo)

        dG = array(dG_nuevo)
        #print "dG despues libre:    ", dG


        dG = self.transformacion_colision(mesa, frontera, t, dG, normas, v_antes, v_despues, theta, alfa)

        #print "dG ya evol:   ", list(dG)

        dG, norms = ortoGS(dG)#; print "dG ortogonalizado:   ", list(dG)
        #print norms
        normas.append(norms)

        return array(dG)





class Orbita():

    """
    Objeto que contiene las trayectroias e informacion de la dinamica
    """

    def __init__(self, mesa):

        """
        Inicializa el objeto orbita, con pos y vel iniciales
        """

        self.mesa = mesa

        self.R = [] #Puntos de colision con la frontera
        self.S = [] #Posicion de colision con la frontera dada en terminos de la long de arco
        self.V = [] #Velocidades
        self.F = [] #Fronteras de colision
        self.T = [0.0] #Tiempo
        self.N = [] #numero de colisiones
        self.t_total = self.T[0]
        self.Angulos_entrada = [] #Angulo de colison
        self.Angulos_salida = [] #Angulo de reflexion
        self.dGamma = []
        self.normas_dGamma = []

        #self.vec_tang_map = rand(2)                #variables para exponentes de lyap con la derivada del mapeo
        #vec_tang_map /= norm(vect_tang_map)

        self.lyap_discrete = []
        self.lyap_flow = []


        self.x0 = array(self.mesa.pos_ini()) #vector de posicion inicial
        self.v0 = array(self.mesa.vel_ini()) #vector de velocida inicial

        self.frontera_previa = -1
        #self.Reglas = ReglasDeColision()
        #self.Mapeos = MapeosTangentes()

        self.alfa = -1 #para saber si no se ha hecho ninguna corrida


        ############ variables que se utilizan para calcular la dimension del atractor en las funciones cubeirta y dimension

        #self.E = -1
        self.centros_de_cubierta = -1
        self.epsilon = -1
        self.M = -1




    def flujo(self, pasos, regla, alfa, calcular_lyapunov=True):

        """
        Continua la dinamica por n colisiones mas utilizando la regla de colision que se pase como funcion al argumento
        "reflexion" dentro de la mesa que se pase como "frontera" (instancia de un objeto). La variable reflex debe ser "a" o "d"
        segun sea mi regla o la de David que se quiera utilizar.
        """

        self.alfa = alfa

        self.Reglas = ReglasDeColision(regla)
        self.Mapeos = MapeosTangentes(regla)

        reflexion = self.Reglas.colision

        #m_l = self.Mapeos.m_l


        ######## Se seleciona la regla de colision para no tener que hacer un if en cada paso

        #if reflex == "d":
            #reflexion = self.Reglas.david
            #Dmap = self.Mapeos.A_d
        #elif reflex == "a":
            #reflexion = self.Reglas.alonso
            #Dmap = self.Mapeos.A_a
        #else:
            #print "Regla de reflexion invalida\nDebe ser: 'a' o 'd'"

        ######## Se selecciona la funcion de pos_de_arco de la mesa dependiendo si es regular o no para evitar un if en cada iteraci贸n
        ######## Esto se podr铆a evitar si se define un clase de poligonos que hereda si es reguar o irregular

        if self.mesa.regular:
            pos_de_arco = self.mesa.pos_de_arco_regular
        else:
            pos_de_arco = self.mesa.pos_de_arco_no_regular


        ############ se inicialza el punto inicial

        x = self.x0
        v = self.v0

        ############ Se inicializan los vectores tangentes Gamma

        #self.dGamma = rand(4,4)
        #self.dGamma, var_que_no_se_usa = ortoGS(self.dGamma)

        self.dGamma = identity(4)

        ########################################################################

        #if len(self.R) == 0: # Si la dinamica ya se habia empezado, la retoma
            #x = self.x0
            #v = self.v0
        #else:
            #x = self.R[-1]
            #v = self.V[-1]

        for i in range(pasos):

            t, self.frontera_previa = self.mesa.tiempo_interseccion(x, v, self)

            self.t_total += t

            x_nuevo = x + v*t

            v_nuevo = reflexion(x, v, self.frontera_previa, self.mesa, alfa)

            self.frontera = self.frontera_previa


            ###### Evolucionamos los vectores tangentes!
            ###### CHECAR QUE EN COLISION SI SEA FRONTERA_PREVIA LA ADECUADA!!!!!!!

            #vuelo_libre = lambda x: dot(m_l(t),x)
            #colision = lambda x: dot(Dmap(self.mesa.normales[self.frontera_previa], v, alfa, self.Reglas.angulo_salida), x)

            #self.dGamma = map(vuelo_libre, self.dGamma)
            #self.dGamma = map(colision, self.dGamma)

            #print "dGamma:\n", self.dGamma

            #self.dGamma = self.Mapeos.evolucion_tangentes(self.mesa, self.frontera, t, self.dGamma, self.normas_dGamma, v, self.Reglas.angulo_salida, self.alfa )

            if calcular_lyapunov:
                self.dGamma = self.Mapeos.evolucion_tangentes(self.mesa, self.frontera, t, self.dGamma, self.normas_dGamma, v, v_nuevo, self.Reglas.angulo_entrada, self.alfa );#print self.dGamma

            #print "angulo:\t", self.Reglas.angulo_salida
            #for i in range(len(self.dGamma)):

                #self.dGamma[i] = dot(Dmap(self.mesa.normales[self.frontera_previa], v, alfa, self.Reglas.angulo_salida), dot(m_l(t), self.dGamma[i]) )


            #self.normas_dGamma.append( map(linalg.norm, self.dGamma) )

            #ortoGS(self.dGamma)

            v = v_nuevo
            x = x_nuevo
            #####################################################################



            #### Actualizamos valores

            self.R.append(x)

            self.V.append(v)

            self.F.append(self.frontera_previa)

            self.T.append(self.t_total)

            self.N.append(i)

            self.Angulos_entrada.append(self.Reglas.angulo_entrada)

            self.Angulos_salida.append(self.Reglas.angulo_salida)

            self.S.append(pos_de_arco(x, self.frontera_previa))



        #self.normas_dGamma = array(self.normas_dGamma)

        self.lyap_discrete = lyapunov(self.normas_dGamma, self.N, 10000)#int(len(self.T)*0.1) + 1 )
        self.lyap_flow = lyapunov(self.normas_dGamma, self.T, 10000)

        #t_lyap_ini = int(len(self.T)*0.3)

        #for i in range(len(self.dGamma)):

            #self.lyap.append( (1.0/(self.T[-1] - self.T[t_lyap_ini]))*sum(log(self.normas_dGamma[t_lyap_ini:][i])) )


    def draw_espacio_fase(self, a_partir_de = 0):

        """
         Dibuja el espacio fase, poniendo puntos a partir de la colison "a_partir_de"
        """

        S = array(self.S)
        Angulos_salida = array(self.Angulos_salida)
        title(u'Espacio Fase (sección de Poincaré)')
        ylabel(u'Ángulo de salida ($\\theta$)')
        xlabel(u'Longitud de arco ($s$)')

        plot(S[a_partir_de:], Angulos_salida[a_partir_de:], ".", markersize = 0.2 )
        axis([0, self.mesa.longitud_de_arco, -pi/2, pi/2])

        show()

    def draw_trayectoria(self):

        """
        Dibuja la trayectoria de la particula (espacio de configuraciones)
        """

        R = array(self.R)
        V = array(self.V)

        title("Trayectoria")

        self.mesa.draw()

        plot(R[:,0], R[:,1], "-om")
        plot(R[0,0], R[0,1], "ok")             #dibuja el punto inicial en negro
        quiver(R[0,0], R[0,1], V[0,0], V[0,1]) #dibuja el vector velocidad en el punto inicial

        show()


    def draw_mapeo(self, a_partir_de = 0):

        """
        Dibuja el "mapeo", solo esta bien definido para el caso slap
        """

        S = array(self.S)
        S1 = S[:-1]
        S2 = S[1:]

        xmax = self.mesa.longitud_de_arco

        X = linspace(0, xmax, 50)

        clf()

        plot(S1[a_partir_de:], S2[a_partir_de:], "o", markersize = 0.9)
        plot(X, X, "b--")
        axis([0, xmax, 0, xmax])

        title("Mapeo")
        xlabel("$s_{n}$")
        ylabel("$s_{n+1}$")

        show()


    def guardar(self, nombre = "orbita.dat"):

        """
        Guarda los datos de la dinamica en un archivo de datos

        Se le pasa como argumento un str con el nombre del archivo
        """

        R = array(self.R)
        S = array(self.S)
        F = array(self.F)
        A = array(self.Angulos_salida)
        A_en = array(self.Angulos_entrada)
        T = array(self.T)

        archivo = open(nombre, 'w')


        header = "#poligono con {0} lados, alfa = {1}, pasos = {2}\n# n\tt\t\tx\ty\t\tf(frontera de colision)\ts(long de arco)\tA_en(angulo de entrada)\tA(angulo de salida)\n\n".format(self.mesa.n, self.alfa, len(R))
        archivo.write(header)


        for i in range(len(self.R)):

            line = "{0:<8g} {1:<13g} {2:<13g} {3:<13g} {4:<4g} {5:<13g} {6:<13g} {7:<13g}\n".format(i, T[i], R[i,0], R[i,1], F[i], S[i], A_en[i], A[i])

            archivo.write(line)

            #print str(i)+"\t"+str(self.R[i][0])+"\t"+str(self.R[i][1])+"\t"+str(self.S[i])+"\t"+str(self.Angulos_salida[i])
        archivo.close()

        print "datos guardados en '{0}' !".format(nombre)


    #def draw_cubierta(self, a_partir_de = 0):

        #"""
        #Grafica la cubierta que se puede utilizar para calcular la dimension fractal del atractor

        #Se tiene que correr antes la funcion dimension
        #"""

        #R = array(self.R).copy()
        #S = array(self.S).copy()
        #A = array(self.Angulos_salida).copy()

        #R = R[a_partir_de:]
        #S = S[a_partir_de:]
        #A = A[a_partir_de:]

        #l = self.centros_de_cubierta
        #diametro = self.epsilon * 2.0

        #Coordenadas = [array([S[i], A[i]]) for i in range(len(R))]


        ## se define el arreglo de los circulos que forman la cubierta
        #E = [ matplotlib.patches.Ellipse(xy = array(Coordenadas[i]), width = diametro, height = diametro) for i in l]

        ##self.E = E

        #fig = figure()
        #ax = fig.add_subplot(111)

        #plt = plot(S, A, ".", markersize = 0.8 )


        #for e in E:
            #ax.add_artist(e)
            #e.set_clip_box(ax.bbox)
            #e.set_alpha(0.3)
            #e.set_facecolor("r")

        #xmax = self.mesa.longitud_de_arco

        #ax.set_xlim(0, xmax)
        #ax.set_ylim(-pi/2, pi/2)



    #def dimension(self, circulos = 1000, epsilon = 0.075, a_partir_de = 0):

        #"""
        #Calcula, la dimension fractal (dimension de informacion) del atractor.

        #"""

        ##E = self.E

        #R = array(self.R).copy()
        #S = array(self.S).copy()
        #A = array(self.Angulos_salida).copy()

        #R = R[a_partir_de:]
        #S = S[a_partir_de:]
        #A = A[a_partir_de:]

        #self.epsilon = epsilon

        #numero_de_pasos = len(R)
        #numero_de_circulos = circulos

        ######################################### Hay que ver como se debe elegir la cubierta para que funcione esto...

        #c = [] # lista con los indices que corresponden al punto que se va a tomar como el centro del circulo....

        #for i in range(numero_de_circulos):
            #num = int( i * ( numero_de_pasos / float(numero_de_circulos ) ))
            #c.append(num)

        #self.centros_de_cubierta = c

        ########################################

        #diametro = 2.0 *epsilon

        #Coordenadas = [array([S[i], A[i]]) for i in range(len(R))]

        ##print Coordenadas


        #P = [] # Esto va a ser la frecuencia relativa con la que un punto cae dentro de la bola, por el momento es un arreglo de ccontadores


        #for cc in c:

            #f = lambda x: x - Coordenadas[cc]
            #contador = 0

            #for x in Coordenadas:

                #if norm(f(x)) < epsilon:
                    #contador += 1

                ##print contador

            ##print P
            #P.append(contador)


        #print P


        #g = lambda x: x / float(len(R))
        #P = map(g, P)

        ##print P

        ##g = lambda x: x*log(x)

        ##h= map(g, P)

        #for i in range(len(P)):
            #if P[i] > 0.0:
                #P[i] = P[i]*log(P[i])
            #else:
                #P[i] = 0



        #H = 0
        #for i in range(len(P)):
            #H += P[i]

        #print "H positivo:\t", H

        #H = -1*H

        #D = H / log(1/epsilon)

        #print "H:\t",H
        #print "log(1/E):\t", log(1/epsilon)
        #print "Dimension de Infromacion:\t", D

        #return D


    #def dim_cap(self, epsilon):

        #"""
        #Calcula la dimension de capacidad (capacity / box-counting dimension)
        #"""

        #### Primero hay que dividir el espacio fase en cajas de "tamano" epsilon

        #S = list(array(self.S))   # Lo hago de esta manera para que al manipular S y A no altere self.S y self.A ...
        #A = list(array(self.Angulos_salida))

        #X = array(column_stack((S,A)))  # este arreglo contiene los vectores de los puntos del espacio fase...
        #print X

        #X_mod = X.copy()

        #s_max = self.mesa.longitud_de_arco
        #s_min = 0

        #a_max = pi / 2
        #a_min = -pi / 2

        #s_rango = s_max - s_min
        #a_rango = a_max - a_min

        #num_cajas = a_rango / epsilon ## este es el numero de cajas que se pueden formar en una "direccion" del espacio fase

        #longitud_de_caja = s_rango / num_cajas
        #altura_de_caja = epsilon

        #print longitud_de_caja, altura_de_caja

        ##centros_s = linspace((longitud_de_caja/2), (s_max - longitud_de_caja/2), num_cajas)
        ##centros_a = linspace(a_min + (altura_de_caja/2), (pi/2 - altura_de_caja/2), num_cajas)

        #centros_s = arange(longitud_de_caja/2.0, s_max, longitud_de_caja)
        #centros_a = arange(-pi/2 + epsilon/2.0, pi/2, epsilon )

        #centros_de_cajas = column_stack((centros_s, centros_a))
        #self.centros_de_cubierta = centros_de_cajas

        #plot(centros_de_cajas[:,0], centros_de_cajas[:,1], "o")

        #contador = 0  # va a ir contando las cajas ocupadas


        #for i in range(len(X)):

            #vectores_posibles = [] # lista de vectores posibles para la primer columna de cajas

            #for j in range(len(centros_s)):   # Me gustaria poder ir quitando los vectores de la lista para hacer mas eficiente el metodo
                #if ((X[i,0] > centros_s[j] - longitud_de_caja) and (X[i,0] < centros_s[j] + longitud_de_caja)):
                    #vectores_posibles.append(X[j])


            #c = 0
            #for k in range(len(centros_a)):
                #for x in list(vectores_posibles):

                    #if ((x[1] > list(centros_a)[k] - altura_de_caja) and (x[1] < list(centros_a)[k] + altura_de_caja)):
                        #c += 1

                #if c > 0:
                    #contador += 1



        ###M = meshgrid(centros_de_cajas[0], centros_de_cajas[1])

        #return contador



    def dimension_hist(self, N):

        """ N_s es el numero de columnas de cajas que se van a tener. Se van a tener el mismo numero de filas de cajas
        en la dimension angular tambien. Asi que el numero de cajas totales van a ser N**2,
        """

        S = array(self.S)   # Lo hago de esta manera para que al manipular S y A no altere self.S y self.A ...
        A = array(self.Angulos_salida)

        A_min = -pi/2.0

        #M = zeros((N,N), dtype='int')
        #self.M = M

        epsilon_s = self.mesa.longitud_de_arco / N
        epsilon_a = pi / N

        #print "epsilon s:\t", epsilon_s
        #print "\nepsilon a:\t", epsilon_a

        s = floor(S / epsilon_s).astype('int') # se dividen las coord entre el tamano de la caja para ver en cual columna/renglon estan...
        a = floor( (A - A_min) / epsilon_a).astype('int')

        cajas = { (s[i], a[i]):0 for i in range(len(s))}   #dict comprehension

        for i in range(len(s)):

            cajas[(s[i], a[i])] += 1

        self.cajas = cajas

        casillas_ocupadas = len(cajas) # Para la dim de capacidad box-counting

        ocupacion = array(cajas.values())
        p = ocupacion/float(sum(ocupacion))

        H = -sum(p*log(p)) #entropia de shannon???
        N_epsilon = casillas_ocupadas #numero de cajas ocupadas

        return epsilon_s, N_epsilon, H



    def itera_dim(self, epsilon_min, guardar = False):


        """
        Itera el procedimiento para calcular la dim de inf, y de capacidad. grafica la curva log-log,
        ajusta una recta y devuelve el resultado

        se le pasa el valor minimo que tomar谩 epsilon
        """

        m_max = int(log2(self.mesa.longitud_de_arco /epsilon_min))

        """  todo este rollo es para que en la grafica log-log se vean equiespaciados
        los puntos en el eje x, esto lo tom茅 del Chua
        """

        rango = []

        for m in range(m_max):
            rango.append(2**m)

        rango = map(int, rango)


        arreglo_epsilon = []
        arreglo_N_epsilon = []
        arreglo_H = []

        for n in rango:

            epsilon, N_epsilon, H = self.dimension_hist(n)

            arreglo_epsilon.append( epsilon )
            arreglo_N_epsilon.append( N_epsilon )
            arreglo_H.append( H )

            print "epsilon = ", epsilon


        arreglo_epsilon = array(arreglo_epsilon)
        arreglo_N_epsilon = array(arreglo_N_epsilon)
        arreglo_H = array(arreglo_H)

        f = lambda x: 1/x

        uno_sobre_epsilon = map(f, arreglo_epsilon)


        self.arreglo_log_uno_sobre_epsilon = log(uno_sobre_epsilon)
        self.arreglo_log_N_epsilon = log(arreglo_N_epsilon)
        self.arreglo_H = arreglo_H

        #clf()

        #title(u"Gráficas para determinar la dimensión fractal del atractor para $ \alpha = {0} $".format(self.alfa))
        #xlabel(r"$log(\frac{1}{\epsilon})$")
        #ylabel(r"$\mathcal{H}_{q}$")

        #grid(True)

        #Dim Capacidad
        #l1 = plot(log(uno_sobre_epsilon), arreglo_log_N_epsilon), ":o", label= u"$\mathcal{H}_0$ (Dimensión de Capacidad)")

        #Dim Informacion
        l2 = plot(log(uno_sobre_epsilon), arreglo_H, ":o", label = u"$\mathcal{H}_1$ (Dimensión de Información)")

        #legend(loc = 4)
        #show()

        if guardar:

            nombre = "alfa_{0}_datos_dim.dat".format(int(100*self.alfa))

            arreglo_dim = column_stack((self.arreglo_log_uno_sobre_epsilon, self.arreglo_log_N_epsilon, self.arreglo_H))

            #arreglo_cap = column_stack((self.arreglo_log_uno_sobre_epsilon, self.arreglo_log_N_epsilon))
            #arreglo_inf = column_stack((self.arreglo_log_uno_sobre_epsilon, self.arreglo_H))

            savetxt(nombre , arreglo_dim) # datos de dim de cap
            #savetxt(nombre + "inf.dat", arreglo_inf ) # datos de dim de inf


#def ortoGS(V):

    #"""
    #Ortogonaliza una lista de vectores en el orden dado, con proceso Gram-Shmidt
    #"""

    #n = len(V)

    #for i in range(len(V)):
        #V[i] /= linalg.norm(V[i])

    #for i in range(1,n):
        #for j in range(0, i):

            #V[i] -= ( dot(V[i], V[j]) / dot(V[j], V[j]) )*V[j]

        #V[i] /= sqrt(dot(V[i], V[i]))


def ortoGS(vector_list):

    """
    Ortogonaliza una lista de vectores en el orden dado, con proceso Gram-Shmidt
    y da como salida una lista de las normas antes de la ultima normalizacion

    NO REGRESA LA LISTA DE VECTORES, SOLO LA ALTERA! Regresa la lista de normas
    """

    #V = sorted(array(Vector_list).copy(), key=norm, reverse=True)
    V = array(vector_list).copy()


    #print "V antes:  ", list(V)

    norms = []
    for j in range(len(V)):

        for i in range(j):
            #print "j,i", j, i

            V[j] -= dot(V[j], V[i])*V[i]
            #print "V[j=", V[j]

        norma_j = norm(V[j])
        V[j] /= norma_j

        norms.append(norma_j)

    #print "V : " , list(V), "\n"
    #importa el orden en el que se ortogonalizan los vectores??  Siempre se debe embezar por el más largo ????????
    #return norms
    return V.copy(), norms

    #norms = []
    #for i in range(len(V)):
        #for j in range(i+1, len(V)):

            #V[j] -= dot(V[i], V[j])*V[i]

        #norma_i = norm(V[i])
        #V[i] /= norma_i

        #norms.append(norma_i)

    #return norms


def lyapunov(normas, tiempos, transientes):

    """
    Toma la lista de las normas de los vectores tangentes y devuevle los exponentes de lyapunov
    asociados con los vectores tangentes.


    El argumento transientes es el numero de puntos que se deben ignorar para asegurar que solo
    se contemplen puntos sobre el atractor.

    Hay que tener mucho cuidado con la forma del arreglo de las normas!!!!
    """

    normas = array(normas)
    lyap = []

    for i in range(len(normas[0])):

        exponente = sum(log(normas[transientes:,i])) / (tiempos[-1] - tiempos[transientes])#(tiempos[-1] - tiempos[transientes])

        lyap.append(exponente)

    return lyap



    #self.normas_dGamma = array(self.normas_dGamma)

    #self.lyap = []

        #t_lyap_ini = int(len(self.T)*0.3)

        #for i in range(len(self.normas_dGamma)):

            #self.lyap.append( (1.0/(self.T[-1] - self.T[t_lyap_ini]))*sum(log(self.normas_dGamma[t_lyap_ini:][i])) )



#p = Poligono([0,0], 3, 1, regular =  True )
#orb = Orbita(p)
#orb.flujo(1000000, "d", 0.65)





#param = arange(1, -0.01, -0.01)
#lyapunov_exponents = []

#for lamb in param:


    #p = Poligono([0,0], 3, 1, regular =  True )
    #orb = Orbita(p)
    #orb.flujo(100000, "d", lamb)

    #lyapunov_exponents.append(orb.lyap)

#lyapunov_exponents = array(lyapunov_exponents)
#plot(param, lyapunov_exponents[:,0])
#plot(param, lyapunov_exponents[:,1])
#plot(param, lyapunov_exponents[:,2])
#plot(param, lyapunov_exponents[:,3])
#show()


#p = Poligono([0,0], 3, 1, regular =  True )
#orb = Orbita(p)
#orb.flujo(1000, "d", 0.001)

#print "Vectores:\n", orb.dGamma
#print "\nNormas:\n", orb.normas_dGamma[-1]
#print "\nExponentes:\n",orb.lyap
#print "\nsuma de exponentes:", sum(orb.lyap)