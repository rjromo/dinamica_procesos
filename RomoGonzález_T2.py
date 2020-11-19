# -*- coding: utf-8 -*-
"""RomoGonzález_t2.ipynb

##RUBÉN ROMO
###TAREA 2
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


opt = int((input("Presione 1 para graficar modelo no lineal \n2 Para graficar modelo lineal \n3 Para graficar mediante solve_ivp")))

if opt ==1:
# MODELO NO LINEAL
  def cstr_NO_lineal (x,t):
      Ca = x[0]
      Cb = x[1]
      Cc = x[2]
      
      qi= 1.1  #m3/min
      K1 = 0.29 #min
      K2 = 0.14  #min
      v = 19 #m3
      Cai=5 #kmol/m3
      Cbi=4.5 #kmol/m3
      
      dCadt = (qi/v)*(Cai-Ca) - (K1*Ca)
      dCbdt = (qi/v)*(Cbi-Cb) + (K1*Ca) - (K2*Cb)
      dCcdt = -(qi/v)*Cc + (K2*Cb)
    
      return [dCadt, dCbdt, dCcdt]

  Cai = 5 # kmol/m3
  Cbi = 4.5 # kmol/m3
  Cci = 0.0 # kmol/m3

  xi = [Cai, Cbi, Cci]

  t0 = 0 #tiempo inicial
  tf = 100 #tiempo final
  puntos = 101 #puntos

  t=np.linspace(t0,tf,puntos)

  y= odeint(cstr_NO_lineal, xi, t)




  plt.title("Ca vs t - no lineal")
  plt.plot(t, y[:,0], 'r')
  plt.xlabel("Time (min)")
  plt.ylabel("Ca (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Cb vs t - no lineal")
  plt.plot(t, y[:,1], 'y')
  plt.xlabel("Time (min)")
  plt.ylabel("Cb (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Cc vs t - no lineal")
  plt.plot(t, y[:,2], 'g')
  plt.xlabel("Time (min)")
  plt.ylabel("Cc (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Ca, Cb, Cc vs t - no lineal")
  plt.plot(t, y[:,0], 'r', label='Ca')
  plt.plot(t, y[:,1], 'y', label='Cb')
  plt.plot(t, y[:,2], 'g', label='Cc')
  plt.xlabel("Time (min)")
  plt.ylabel("C*(kmol/m3)")
  plt.legend(loc='best')
  plt.grid()
  plt.show()

if opt ==2:
  # MODELO LINEAL
  def cstrlineal(x,t):
      ca_desv = x[0]
      cb_desv = x[1]
      cc_desv = x[2]
      
      qi= 1.1  #m3/min
      K1 = 0.29 #min
      K2 = 0.14  #min
      v = 19 #m3
      Cai=5 #kmol/m3
      Cbi=4.5 #kmol/m3
      Cc=0 #kmol/m3

      # steady state
      Ca0= 0.83  #kmol/m3  
      Cb0= 2.58  #kmol/m3
      Cc0= 6.09 #kmol/m3    
      qi0 = qi
      
      
      #expresiones derivadas
      k11 = (-qi/v)-K1
      k12 = 0
      k13 = 0
      b1 = (Cai - Ca0)/v
      k21 = K1
      k22 = (-qi/v)-K2
      k23 = 0
      b2 = (Cbi - Cb0)/v
      k31 = 0
      k32 = K2
      k33 = -(qi/v)
      b3 = (-Cc0)/v

    #EDOs
      dCadt = k11*(ca_desv-Ca0) + k12*(cb_desv-Cb0) + k13*(cc_desv-Cc0) + b1*(qi-qi0)
      dCbdt = k21*(ca_desv-Ca0) + k22*(cb_desv-Cb0) + k23*(cc_desv-Cc0)+ b2*(qi-qi0)
      dCcdt = k31*(ca_desv-Ca0) + k32*(cb_desv-Cb0) + k33*(cc_desv-Cc0) + b3*(qi-qi0)
      
      return [dCadt,dCbdt,dCcdt]

  t0 = 0 #tiempo inicial
  tf = 100 #tiempo final
  puntos = 101 #puntos

  t2=np.linspace(t0,tf,puntos)  #vector de tiempo

  Cai2=5 #kmol/m3
  Cbi2=4.5 #kmol/m3
  Cci2=0 #kmol/m3

  xi2=[Cai2,Cbi2,Cci2]  #vector valores iniciales

  y2=odeint(cstrlineal, xi2, t2) #resolviendo ODEs


  plt.title("Ca vs t - lineal")
  plt.plot(t2, y2[:,0], 'r--')
  plt.xlabel("Time (min)")
  plt.ylabel("Ca (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Cb vs t - lineal")
  plt.plot(t2, y2[:,1], 'y--')
  plt.xlabel("Time (min)")
  plt.ylabel("Cb (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Cc vs t - lineal")
  plt.plot(t2, y2[:,2], 'g--')
  plt.xlabel("Time (min)")
  plt.ylabel("Cc (kmol/m3)")
  plt.grid()
  plt.show()

  plt.title("Ca,Cb,Cc vs t - lineal")
  plt.plot(t2, y2[:,0], 'r', label='Ca')
  plt.plot(t2, y2[:,1], 'y', label='Cb')
  plt.plot(t2, y2[:,2], 'g', label='Cc')
  plt.xlabel("Time (min)")
  plt.ylabel("C*(kmol/m3)")
  plt.legend(loc='lower right')
  plt.grid()
  plt.show()




from scipy.integrate import solve_ivp

if opt ==3:

  # MODELO NO LINEAL SOLVE_IVP
  def cstr_NO_LINEAL_IVP (t,x):
      Ca = x[0]
      Cb = x[1]
      Cc = x[2]
      
      qi= 1.1  #m3/min
      K1 = 0.29 #min
      K2 = 0.14  #min
      v = 19 #m3
      Cai=5 #kmol/m3
      Cbi=4.5 #kmol/m3
      
      dCadt = (qi/v)*(Cai-Ca) - (K1*Ca)
      dCbdt = (qi/v)*(Cbi-Cb) + (K1*Ca) - (K2*Cb)
      dCcdt = -(qi/v)*Cc + (K2*Cb)
    
      return [dCadt, dCbdt, dCcdt]

  Cai = 5 # kmol/m3
  Cbi = 4.5 # kmol/m3
  Cci = 0.0 # kmol/m3
  xi = [Cai, Cbi, Cci]

  ts = [0,100]
  t0 = 0 #tiempo inicial
  tf = 100 #tiempo final
  puntos = 101 #puntos
  t = np.linspace(t0,tf,puntos)

  y= solve_ivp(cstr_NO_LINEAL_IVP, ts, xi, method="RK45", t_eval=t)
  yt = y.y.T


  plt.title("Modelo no lineal")
  plt.plot(y.t, yt[:,0], 'r', label='Ca')
  plt.plot(y.t, yt[:,1], 'y', label='Cb')
  plt.plot(y.t, yt[:,2], 'g', label='Cc')
  plt.xlabel("Time (min)")
  plt.ylabel("C*(kmol/m3)")
  plt.legend(loc='best')
  plt.grid()
  plt.show()
  

  