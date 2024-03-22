import numpy as np
import matplotlib.pyplot as plt
import time

# ----- Inputs ----- #

"""
input_vel = input("Inserire un valore input per selezionare la velocità:\n "
			    "0: c(x) = 1.3\n "
				"1: c(x) = -0.7\n "
				"2: c(x) = x\n "
				"3: c(x) = x^2\n "
				"4: c(x) = atan(x)\n "
				"5: c(x) = sin(x)\n")
input_u0 = input("Inserire un valore input per selezionare il dato iniziale:\n "
			    "0: u0(x) = sin(2πx)\n "
				"1: u0(x) = sin(4πx)\n "
				"2: u0(x) = 0\n "
				"3: u0(x) = x^2\n "
				"4: u0(x) = 1-x^2\n "
				"5: u0(x) = 1-|x|\n")
xmin = input("Inserire un valore input float per xmin: ")
xmax = input("Inserire un valore input float per xmax: ")
T = input("Inserire un valore input float per T: ")
Nx = input("Inserire un valore input float per Nx (numero di passi nello array x): ")
input_cond = input("Inserire un valore input per selezionare le condizioni sul bordo:\n "
				   "0: Condizione periodiche\n "
				   "1: Condizione al bordo inflow g(t)=sin(t)\n "
				   "2: Condizione Direchlet omogenea sul bordo inflow\n "
"""

# ----- Funzioni e parametri della griglia ----- #

xmin, xmax, T = -1, 3, 3
Nx = int(5e2)
x = np.linspace(xmin, xmax, Nx)
dx = (xmax-xmin)/(Nx-1)
print("Valore di dx: ",dx)

def c(x, input):
	if input == 0:
		return 1.3
	elif input == 1:
		return -0.7
	elif input == 2:
		return x
	elif input == 3:
		return x**2
	elif input == 4:
		return np.arctan(x)
	elif input == 5:
		return np.sin(x)

c = c(x, 1)
dt = dx/abs(c)
lam = dt/dx
Nt = int(T/dt + 1)

def u0(x, input, Nx = Nx):
	u0 = np.zeros(Nx)
	if input == 0:
		mask = (-1<=x)*(x<=1)
		u0 = np.sin(2*np.pi*x)*mask
	elif input == 1:
		mask = (-1<=x)*(x<=0)
		u0 = np.sin(4*np.pi*x)*mask
	return u0


def UPWIND(Um, c, n, Nx = Nx, Nt = Nt, lam = lam):
    Up = Um[:]-c*lam*(Um[R]-Um[L])/2+abs(c)*lam*(Um[R]-2*Um[:]+Um[L])/2
    if n==Nt-2:
    	msg = "UPWIND"
    	return Up, msg
    return Up, None

def LAX_FRIEDRICHS(Um, c, n, Nx = Nx, Nt = Nt, lam = lam):
	Up = (Um[R]+Um[L])/2-c*lam/2*(Um[R]-Um[L])
	if n==Nt-2:
		msg = "LAX FRIEDRICHS"
		return Up, msg
	return Up, None

def LAX_WENDROFF(Um, c, n, Nx = Nx, Nt = Nt, lam = lam):
	Up = Um[:]-c*lam/2*(Um[R]-Um[L])+c**2/2*lam**2*(Um[R]-2*Um[:]+Um[L])
	if n==Nt-2:
		msg = "LAX WENDROFF"
		return Up, msg
	return Up, None

def condp(a, b, z, c):
    l=b-a
    if c>0:
    	maskm = (z<a)
    	z += l*maskm
    if c<0:
    	maskp = (z>b)
    	z -= l*maskp
    return z

I = np.arange(Nx)
L = I - 1
R = I + 1
# Condizioni periodiche
L[0] = Nx-1
R[-1] = 0
# Condizioni al bordo inflow


U = np.zeros((Nx,Nt)) # Soluzione numerica
Uex = np.zeros((Nx,Nt)) # Soluzione essata
U[:,0] = u0(x, 0)
Uex[:,0] = u0(x, 0)
# Erorri
err_l1 = np.zeros(Nt-1)
err_linf = np.zeros(Nt-1)

# ----- Ciclo principale ----- #

start_time = time.time()
for n in range(Nt-1):
    U[:,n+1], msg = LAX_WENDROFF(U[:,n], c, n)
    Uex[:,n+1] = u0(condp(xmin,xmax,x-c*n*dt,c),0)
    U_l1 = np.linalg.norm(U[:,n+1],1)
    U_linf = np.linalg.norm(U[:,n+1],np.inf)
    Uex_l1 = np.linalg.norm(Uex[:,n+1],1)
    Uex_linf = np.linalg.norm(Uex[:,n+1],np.inf)
    err_l1[n] = abs((U_l1-Uex_l1)/Uex_l1)
    err_linf[n] = abs((U_linf-Uex_linf)/Uex_linf)
end_time = time.time()
print("Tempo di computazione della U numerica e degli errori: ", end_time-start_time, "s")
print("Errori metodo "+str(msg)+":\n"
	  "	Errori norma l1: "+str(np.mean(err_l1))+"\n"
	  "	Errori norma linf: "+str(np.mean(err_linf)))

# ----- Plot ----- #

plt.figure()
line, = plt.plot(x, U[:,0], "b.", label="Soluzione Numerica ("+str(msg)+")")
line2, = plt.plot(x, Uex[:,0], "r-", label="Soluzione Essata")
plt.legend()
plt.title("Problemi iperbolici (Trasporto)")
plt.xlim(xmin, xmax)
plt.xlabel("x")
plt.ylabel("U")
#plt.plot(x, u0(x, 0), "k-", label="U0(x)")
plt.legend()
for n in range(Nt - 1):
    line.set_ydata(U[:,n+1])
    line2.set_ydata(Uex[:,n+1])
    plt.pause(0.001)
plt.show()