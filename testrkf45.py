#Parham Adiban 1000639446

'''
Using a first degree ODE as a test case in order test the RKF45 algorithm.
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

# Analytic solution to the ODE
def d(t):
	return exp(-2*t)

s = [r]
analytic = [d(a)]


b21 =   2.500000000000000e-01  #  1/4
b31 =   9.375000000000000e-02  #  3/32
b32 =   2.812500000000000e-01  #  9/32
b41 =   8.793809740555303e-01  #  1932/2197
b42 =  -3.277196176604461e+00  # -7200/2197
b43 =   3.320892125625853e+00  #  7296/2197
b51 =   2.032407407407407e+00  #  439/216
b52 =  -8.000000000000000e+00  # -8
b53 =   7.173489278752436e+00  #  3680/513
b54 =  -2.058966861598441e-01  # -845/4104
b61 =  -2.962962962962963e-01  # -8/27
b62 =   2.000000000000000e+00  #  2
b63 =  -1.381676413255361e+00  # -3544/2565
b64 =   4.529727095516569e-01  #  1859/4104
b65 =  -2.750000000000000e-01  # -11/40

# Coefficients used to compute local truncation error estimate.  These
# come from subtracting a 4th order RK estimate from a 5th order RK
# estimate.

r1  =   2.777777777777778e-03  #  1/360
r3  =  -2.994152046783626e-02  # -128/4275
r4  =  -2.919989367357789e-02  # -2197/75240
r5  =   2.000000000000000e-02  #  1/50
r6  =   3.636363636363636e-02  #  2/55

# Coefficients used to compute 4th order RK estimate

c1  =   1.157407407407407e-01  #  25/216
c3  =   5.489278752436647e-01  #  1408/2565
c4  =   5.353313840155945e-01  #  2197/4104
c5  =  -2.000000000000000e-01  # -1/5

# Set t and x according to initial condition and assume that h starts
# with a value that is as large as possible.

t = a
tol = 1e-5

#Implement the runge-kutta-fehlberg method
while t<b:

	if t + h > b:
		h = b - t;

	k1 = h * f(r)
	k2 = h * f(r + b21 * k1)
	k3 = h * f(r + b31 * k1 + b32 * k2)
	k4 = h * f(r + b41 * k1 + b42 * k2 + b43 * k3)
	k5 = h * f(r + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
	k6 = h * f(r + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)

	z = abs((r1*k1 + r3*k3 + r4*k4 + r5*k5 + r6*k6)/h)

	if z <= tol:
		t += h
		r += c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
		s.append(r)
		analytic.append(d(t))


	h = h * min(max((tol/2/z)**0.25, .1), 4.0)

# Plot the Numerical method vs Analytical method
figure()
title('Test of RKF45 method')
plot(s, 'o', label = 'Numerical')
plot(analytic, label = 'Analytic')
legend()
show()

error = (abs(array(analytic) - array(s)))/array(analytic)*100

print(error)