#Parham Adiban 1000639446
'''
Using a first degree ODE as a test case in order test the Adams-Moulton algorithm.
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
s = []
analytic = []

# Coefficient of Adam-Moulton method
b1 = 3/8
b2 = 19/24
b3 = -5/24
b4 = 1/24

c1 = 55/24
c2 = -59/24
c3 = 37/24
c4 = -9/24

# Find analytical solution
for t in tpoints:
	analytic.append(d(t))

# Apply runge kutta to find first 3 values
for i in range(3):
	s.append(r)
	r = r + rk4(r)

# Apply Adam-Moulton method
for j in range(3, N):
	s.append(r)

	# Weight function
	w = r + h*(c1*f(s[j]) + c2*f(s[j-1]) + c3*f(s[j-2]) + c4*f(s[j-3]))
	# Position calculated using weight function
	r = r + h*(b1*f(w) + b2*f(s[j]) + b3*f(s[j-1]) + b4*f(s[j-2]))

# Plot the Numerical method vs Analytical method
figure()
title('Test of Moulton method')
plot(s, 'o', label = 'Numerical')
plot(analytic, label = 'Analytic')
legend()
show()

# Calculate the error 
error = (abs(array(analytic) - array(s)))/array(analytic)*100

print(error)
