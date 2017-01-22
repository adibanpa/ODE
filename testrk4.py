#Parham Adiban 1000639446

'''
Using a first degree ODE as a test case in order test the RK4 algorithm.
'''

from __future__ import division, print_function
from pylab import *

# Define constants
a = 0 # Start time
b = 10 # End time
N = 50 # Number of samples
h = (b-a)/N # Timestep

#Time array
tpoints = arange(a, b, h)

# ODE that needs to be solved dy/dt = -2y
def f(y):
	return -2*y

# Initial condition of y(0) = 1
r = 1

# Function to apply the runge-kutta method
def rk4(r_, h=h):
	k1 = h*f(r_)
	k2 = h*f(r_+0.5*k1)
	k3 = h*f(r_+0.5*k2)
	k4 = h*f(r_+k3)
	return (k1+2*k2+2*k3+k4)/6

# Analytic solution to the ODE
def d(t):
	return exp(-2*t)

# Store numerical and analytical solution
numerical = []
analytic = []

# Perform calculations for all time in the time array
for t in tpoints:
	analytic.append(d(t))
	numerical.append(r)
	r += rk4(r)

# Plot the Numerical method vs Analytical method
figure()
title('Test of RK4 method')
plot(numerical, 'o', label = 'Numerical')
plot(analytic, label = 'Analytic')
legend()
show()

# Calculate the error 
error = (abs(array(analytic) - array(numerical)))/array(analytic)*100

print(error)