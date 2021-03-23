# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer, 
# et for hvert intervall mellom to festepunkter, slik at baade banen y, 
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert) 
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal 
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

#Importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

# CONSTANTS
g = 9.81
c = 2/5
m = 0.031

#Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])
#Skriv inn y-verdiene til banens 8 festepunkter i tabellen yfast.


yfast = np.genfromtxt("./data/y_values.txt", delimiter=";", skip_header=2, usecols=1)
experimental_data_x = np.genfromtxt("./data/3_data.txt", delimiter=";", skip_header=2, usecols=1)
experimental_data_v = np.genfromtxt("./data/3_data.txt", delimiter=";", skip_header=2, usecols=5)


#Et vilkaarlig eksempel:

#Erstatt med egne tallverdier avlest i tracker.
#Programmet beregner de 7 tredjegradspolynomene, et
#for hvert intervall mellom to festepunkter,
#med funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')
#Funksjonen cs kan naa brukes til aa regne ut y(x), y'(x) og y''(x)
#for en vilkaarlig horisontal posisjon x, eventuelt for mange horisontale
#posisjoner lagret i en tabell:
#cs(x)   tilsvarer y(x)
#cs(x,1) tilsvarer y'(x)
#cs(x,2) tilsvarer y''(x)
#Her lager vi en tabell med x-verdier mellom 0 og 1.4 m
xmin = 0.000
xmax = 1.401
dx = 0.001
x = np.arange(xmin, xmax, dx)   
#Funksjonen arange returnerer verdier paa det "halvaapne" intervallet
#[xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
#dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400, 
#dvs x blir en tabell med 1401 elementer
Nx = len(x)
y = cs(x)       #y=tabell med 1401 verdier for y(x)
dy = cs(x,1)    #dy=tabell med 1401 verdier for y'(x)
d2y = cs(x,2)   #d2y=tabell med 1401 verdier for y''(x)


def vy(y, c, g):
  new_list = []
  y0 = y[0]

  for elem in y:
    vy_t = np.sqrt((2*g*(y0 - elem)) / (1 + c))
    new_list.append(vy_t)
  return new_list


def k(dy, d2y):
  new_list = []
  
  for i in range(len(dy)):
    k_t = d2y[i] / (1 + (dy[i]**2))**(3/2)
    new_list.append(k_t)
  return new_list


def a(v, k):
  v_squared = [x**2 for x in v]
  v_np = np.array(v_squared)
  k_np = np.array(k)
  new_array = np.multiply(v_np, k_np).tolist()
  return new_array

def beta(dy):
  new_list = []
  for elem in dy:
    new_list.append(np.arctan(elem))
  return new_list

def normal_force(m, g, beta, centr_f):
  new_list = []
  for i in range(len(beta)):
    new_list.append(m*(g*np.cos(beta[i]) + centr_f[i]))
  return new_list

def f(c, m, g, beta):
  new_list = []
  for i in range(len(beta)):
    element = (c * m * g * np.sin(beta[i]))/(1 + c)
    new_list.append(element)
  return new_list

def rel_friction_normal_f(f, normal_f):
  new_list = []
  for i in range(len(f)):
    new_list.append(np.abs(f[i] / normal_f[i]))
  return new_list

def delta_t(v, delta_x, beta):
  v_step = []
  tail = v[0] * np.cos(beta[0])

  for i in range(1, len(v)):
    v_x_n = v[i] * np.cos(beta[i])
    v_step.append(0.5*(tail + v_x_n))
    tail = v_x_n
  
  new_list = []
  for elem in v_step:
    new_list.append(delta_x / elem)
  return new_list

def t_n(time_step):
  new_list = []
  tail = time_step[0]
  new_list.append(tail)

  for i in range(1, len(time_step)):
    new_elem = tail + time_step[i]
    new_list.append(new_elem)
    tail = new_elem
  return new_list

def velocity_dt(dt, dx):
  new_list = []
  for elem in dt:
    new_list.append(dx/elem)
  return new_list

velocity = vy(y, c, g)
curvature = k(dy, d2y)
centripetal = a(velocity, curvature)
inclination = beta(dy)
normal_f = normal_force(m, g, inclination, centripetal)
friction = f(c, m, g, inclination)
rel_n_f = rel_friction_normal_f(friction, normal_f)
t_step = delta_t(velocity, dx, inclination) # t_step shouldn't be plotted
t = t_n(t_step)
v_dt = velocity_dt(t_step, dx)


#Plotteeksempel: Banens form y(x)
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x, rel_n_f)
# plt.plot(x, velocity)
# plt.plot(experimental_data_x, experimental_data_v, color='orange')

# plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$|f/N|$',fontsize=20)
plt.grid()
plt.show()
#Figurer kan lagres i det formatet du foretrekker:
#baneform.savefig("baneform.pdf", bbox_inches='tight')
#baneform.savefig("baneform.png", bbox_inches='tight')
#baneform.savefig("baneform.eps", bbox_inches='tight')