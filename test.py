import numpy as np
import matplotlib.pyplot as plt
import time


#INPUTS:
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

xmin, xmax, T = 0, 7, 5
Nx = int(1e3)
x = np.linspace(xmin, xmax, Nx)
dx = (xmax-xmin)/(Nx-1)

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

def UPWIND(U, c, Nx = Nx, Nt = Nt, lam = lam):
	for n in range(Nt-1):
		#U[:,n+1] = U[:,n]-c*lam*(U[:,n]-U[L,n])
		U[:,n+1] = U[:,n]-c*lam*(U[R,n]-U[L,n])/2+abs(c)*lam*(U[R,n]-2*U[:,n]+U[L,n])/2
	return U
def LAX_FRIEDRICHS(U, c, Nx = Nx, Nt = Nt, lam = lam):
	for n in range(Nt-1):
		#U[:,n+1] = U[:,n]-c*lam*(U[:,n]-U[L,n])
		U[:,n+1] = U[:,n]-c*lam*(U[R,n]-U[L,n])/2+abs(c)*lam*(U[R,n]-2*U[:,n]+U[L,n])/2
	return U
def LAX_WENDROFF(U, c, Nx = Nx, Nt = Nt, lam = lam):
	for n in range(Nt-1):
		#U[:,n+1] = U[:,n]-c*lam*(U[:,n]-U[L,n])
		U[:,n+1] = U[:,n]-c*lam*(U[R,n]-U[L,n])/2+abs(c)*lam*(U[R,n]-2*U[:,n]+U[L,n])/2
	return U


"""
def LAX_FRIEDRICHS():

def LAX_WENDROFF():
"""

I = np.arange(Nx)
L = I - 1
R = I + 1
# Condizioni periodiche
#L[0] = Nx-1
#R[-1] = 0
# COndizioni al bordo inflow



U = np.zeros((Nx,Nt))
start_time = time.time()
U[:,0] = u0(x, 0)
U0 = U[:,0]
end_time = time.time()
print("Computing time for u0 calculation", end_time-start_time)
print(U[:,0])
U = UPWIND(U, c)
print(U)


plt.figure()
line, = plt.plot(x, U0, "b-", label="Soluzione Numerica (UPWIND)")
plt.legend()
plt.title("Problemi iperbolici (Trasporto)")
plt.xlim(xmin, xmax)
plt.xlabel("x")
plt.ylabel("U")
#plt.plot(x, u0(x, 0), "k-", label="U0(x)")
plt.legend()
for n in range(Nt - 1):
    line.set_ydata(U[:,n+1])
    plt.pause(0.001)
plt.show()