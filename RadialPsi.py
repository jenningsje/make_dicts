import numpy as np
import matplotlib.pyplot as plt
from scipy.special import genlaguerre
from scipy.integrate import simps
from scipy.misc import derivative

a = .529
numValues = 20000

def R(n, l, r):
    '''Returns the radial wavefunction for a specific value of r where n is the principle quantum number and l is the angular momentum quantum number.'''
    # Generates an associated Laguerre polynomial as a scipy.special.orthogonal.orthopoly1d
    Lag = genlaguerre(n - l - 1, 2 * l + 1)
    # Initialize the wavefunction value at r to zero
    y = 0.0
    # Calculate the radial wavefunction for the given r
    for i in range(len(Lag) + 1):
        i = float(i)
        # The main equation (What it's all about)
        y += (((r / (n * a)) ** (l)) * (np.exp(-r / (a * n))) * Lag[int(i)] * (((2 * r) / (a * n)) ** i))
    return y



def K(n, l, r):
    '''Returns the radial wavefunction for each value of r where n is the principle quantum number and l is the angular momentum quantum number.'''
    y = np.zeros(len(r))  # Create an array to store the results
    
    # Generates an associated Laguerre polynomial as a scipy.special.orthogonal.orthopoly1d
    Lag = genlaguerre(n - l - 1, 2 * l + 1)
    
    for i in range(len(Lag)):
        # Define a function for R(n, l, r) to use with derivative
        def K_func(x):
            return (((x / (n * a))**(l)) * (np.exp(-x / (a * n))) * Lag[i] * (((2 * x) / (a * n))**i))
        
        # Compute the second derivative of R_func with respect to x
        y += [derivative(K_func, r_val, dx=1e-6, n=2) for r_val in r]
    
    return y

#Approximates integral value and divides by the absolute value.
def normalise(x,y):#Doesnt normalise. Calculates integral, we want area under the graphs. Haven't looked at this in a while, needs checking.
	integral = simps(np.absolute(y),x)
	print(integral, "simps")
	print(type(integral))
	y = y/np.absolute(integral)
	return y

def plotPsi(r,psi):
	plt.subplot(311)
	plt.plot(r, psi)
	plt.grid(True)
	plt.xlabel("r/a (m)")
	plt.ylabel("Psi")

def plotPsiSquared(r,psi):
	psiSquared = psi**2
	psiSquared = normalise(r,psiSquared)
	
	plt.plot(r, psiSquared)
	plt.grid(True)
	plt.xlabel("r/a (m)")
	plt.ylabel("Psi^2")

def plotRadialDistribution(r,psi):
	radialDistribution = 4*np.pi*(r**2)*psi**2
	radialDistribution = normalise(r, radialDistribution)

	plt.plot(r, radialDistribution)
	plt.grid(True)
	plt.xlabel("r/a (m)")
	plt.ylabel("4*pi*(r^2)*Psi")

def plotAll(r,psi):
	plt.figure(1)
	plt.subplot(311)
	plotPsi(r,psi)
	plt.subplot(312)
	plotPsiSquared(r,psi)
	plt.subplot(313)
	plotRadialDistribution(r,psi)

def graphs(r, psi,choice):
	if choice == 1:
		plotPsi(r,psi)
	elif choice == 2:
		plotPsiSquared(r,psi)
	elif choice == 3:
		plotRadialDistribution(r,psi)
	else:
		plotAll(r,psi)

def main():
	numPsi = int(input("How many wavefunctions do you want to draw?"))
	#Note that the whole wavefunction should be drawn or the normalise function will not work correctly
	width = float(input("To what value of r (in units of a) do you want to plot?")) * a
	print("1: Plot the wavefunction")
	print("2: Plot the wavefunction squared (the probability density)")
	print("3: Plot the wavefunction squared times pi r^2 (the radial distribution function)")
	print("4: Plot all 3")
	choice = int(input("Choose an option:"))
	for i in range(numPsi):
		print("Wavefunction " + str(i+1))
		n = float(input("Enter n:"))
		l = float(input("Enter l:"))
		r = np.linspace(0,width, numValues)
		psi = R(n,l,r)

		#Converts r into units of a
		r = r/a
		psi = normalise(r, psi)
		
		graphs(r, psi,choice)
	plt.show()


if __name__ == '__main__':
    main()