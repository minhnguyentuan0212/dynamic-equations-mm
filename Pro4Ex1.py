import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import math
#IMPLICIT METHOD SOLUTION

#configuration
R = []
J = []
t = []
n = 40
#initial value
R0 = 0
J0 = 0
h = 0.1
t0 = 0

#iteration
for i in range(0, n):
  t1 = t0 + h
  R1 = (9*R0 + J0 - 7*h*t1)/8.2
  J1 = -(R0 - 9*J0 + -19*h*t1)/8.2
  R.append(R1)
  J.append(J1)
  t.append(t0)
  R0 = R1
  J0 = J1
  t0 = t1

#plot solution
plt.plot(t, R, 'b', label = "R(t)")
plt.plot(t, J, 'g', label = "J(t)")
plt.xlabel("t")
plt.legend()
plt.legend(loc='best')
plt.grid()
plt.title("Implicit Euler solution")
plt.show()


#EXACT SOLUTION
y0 = [0, 0] #vector [R, J]
t = np.linspace(0, 4, 61)
def pend(y, t,):
    R, J = y
    dydt = [R + J - t, -R + J + 2*t]
    return dydt



#plot solution
sol = odeint(pend, y0, t)
plt.plot(t, sol[:, 0], 'b', label='R(t)')
plt.plot(t, sol[:, 1], 'g', label='J(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.title("Exact Solution (with help of scipy library in Python)")
plt.show()