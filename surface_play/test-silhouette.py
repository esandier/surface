# 

from silhouette import Surface
import time

T0 = time.perf_counter()

#surf = Surface('u','v', 'log(cos(v)/cos(u))', 'u v')
#surf.triangulate(bounds = (-1.5,1.5,-1.5,1.5), res = 53)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(1.2+ cos(u))*cos(v)','(1.2 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (0,6,0, 5.5), res = 147, quotient=('no','no'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

#surf = Surface(
#"u*cos(v)", "u*sin(v)", "v", "u v"
#)  # demande une perturbation pas trop grande à azimuth 0.
#surf.triangulate(bounds=(-5, 5, 0, 25), res=200)
#surf.set_axis(elev=10, azim=8)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u', 'v', 'u*u*u + v*u', 'u v')
#surf.triangulate(bounds = (-1,1,-1,1), res = 32)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*v', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 61)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','2*sin(u*v)', 'u v' ) # fonctionne avec res = 600, res ligne = 200, bornes -5,5,-5,5
#surf.triangulate(bounds = (-5,5,-5,5), res = 700)
#surf.set_axis(elev = 30, azim = 38)
#surf.traitement()
#surf.plot(T0)
## mauvais résultat pour certains points de vue

#surf = Surface('u','v','3*sin(u)*sin(v)', 'u v' ) # marche bien à 200
#surf.triangulate(bounds = (0,15,0,15), res = 203)
#surf.set_axis(elev = 9, azim = 67)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','(1 - (u*u + v*v))**2', 'u v' )  # marche à 100, elev 20, azim 28
#surf.triangulate(bounds = (-1.3,1.3,-1.3,1.3), res = 33)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','3*(sin(u)+sin(v))**2', 'u v' )
#surf.triangulate(bounds = (0,20,0, 20), res = 506) # à (0,20,0,20) res 500, elev 25, azim 35 ça marche pas, même à 1500
#surf.set_axis(elev = 25, azim = 35)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*u', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 51)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','u*u - v*v', 'u v' )
#surf.triangulate(bounds = (-2,2,-2,2), res = 29)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

# surf = Surface('u','v','4*sin(u)', 'u v', bounds = (0,10,-2,2))
# surf.triangulate(res = 43)
# surf.set_axis([1,1,0],[-1,1,3])
# surf.traitement()
# surf.plot_for_terminal(T0)

#surf = Surface('u*(4+ sin(10*exp(-5*(u*u+v*v))))/sqrt(u*u+v*v+.0001)','3*(u*u+v*v)','v*(4+ sin(10*exp(-5*(u*u+v*v))))/sqrt(u*u+v*v+.0001)', 'u v' )
#surf.triangulate(bounds = (-1,1,-1,1), res = 367)
#surf.set_axis(elev = 20, azim = 28)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('u','v','2*sin(u*v)', 'u v' )
#surf.triangulate(bounds = (-4,4,-4,4), res = 176)
#surf.set_axis(elev = 31, azim = 39)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(3+ u*cos(3*v))*cos(v)','(3+ u*cos(3*v))*sin(v) ','u*sin(3*v)', 'u v')
#surf.triangulate(bounds = (-1,1,0, 6.28318530717958647), res = 143, quotient=('no','cy'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

# surf = Surface('u','v','sin(u*u + v*v)', 'u v', bounds = (-5,5,-5,5)) # temps à battre : entre 10 et 12,5 secondes
# surf.triangulate(res = 300) # 600 OK avec x*x/(1+x), même pour -6,6,-6,6. 800 pour x^{3/2}/(1+sqrt(x))
# surf.set_axis([-1,1,0],[-1,-1,2])
# surf.traitement()
# surf.plot_for_terminal(T0)

# surf = Surface('cos(u)*cosh(v)','sin(u)*cosh(v)','v', 'u v', bounds = (0,6.28318,-1.3,1.3), quotient=('cy','no')) # ça déconne complètement
# surf.triangulate(res = 300) 
# surf.set_axis([-1,1,0],[-1,-1,2])
# surf.traitement()
# surf.plot_for_terminal(T0)


#################################################################################################################
#################################################################################################################
#################################################################################################################


#surf = Surface('(3+ cos(u))*cos(v)','(3 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (.3,6, 0,6.2831853071795864769252), res = 100, quotient=('no','cy'))
#surf.set_axis(elev = 30, azim = -30)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(2.5+ cos(u))*cos(v)','(2.5 + cos(u))*sin(v) ','sin(u)', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 0,6.2831853071795864769252), res = 200, quotient=('cy','cy'))
#surf.set_axis(elev = 26.6, azim = 14.88)
#surf.traitement()
#surf.plot(T0)

#surf = Surface('(2+ sin(v))*cos(u)','(2+sin(v))*sin(u) ','v', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 1,10), res = 200, quotient=('cy','no'))
#surf.set_axis(elev = 10, azim = 40)
#surf.traitement()
#surf.plot(T0)

surf = Surface('(1.2+ sin(2*v))*cos(u)','(1.2+sin(2*v))*sin(u) ','v', 'u v', bounds = (0,6.2831853071795864769252, 0,10), quotient=('cy','no')) 
surf.triangulate(res = 300)
surf.set_axis([-1,1,0],[-1,-1,2])
surf.traitement()
surf.plot_for_terminal(T0)

#surf = Surface('cos(u)','sin(u) ','v', 'u v')
#surf.triangulate(bounds = (0,6.2831853071795864769252, 3,9), res = 200, quotient=('cy','no'))
#surf.set_axis(elev = 10, azim = 40)
#surf.traitement()
#surf.plot(T0)

# surf = Surface("(2+ u*cos(v/2))*cos(v)", "(2+ u*cos(v/2))*sin(v) ", "u*sin(v/2)", "v u", bounds=(0, 6.28318, -1, 1), quotient=("mo", "no")) # deconne elev=40, azim=40 si res = 600, uniquement si raccord "mo, no"
# surf.triangulate(res=100)
# surf.set_axis([-1,1,0],[-1,-1,2])
# surf.traitement()
# surf.plot_for_terminal(T0)

#surf = Surface('(2+ u*cos(3*v/2))*cos(v)','(2+ u*cos(3*v/2))*sin(v) ','u*sin(3*v/2)', 'u v')
#surf.triangulate(bounds = (-1.4,1.4,0, 6.2831), res = 100, quotient=('no','mo'))
#surf.set_axis(elev = 8, azim = 50)
#surf.traitement()
#surf.plot(T0)
