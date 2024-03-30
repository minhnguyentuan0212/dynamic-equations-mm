import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import math 
#IMPLICIT METHOD SOLUTION
#configuration
R = [2]
J = [4]
t = [0]
n = 100
#initial value
R0 = 2
J0 = 4
h = 0.1
t0 = 0

#iteration
def quadric(a, b, c):
  dis = b * b - 4 * a * c 
  sqrt_val = math.sqrt(abs(dis)) 
  x1 = (-b + sqrt_val)/(2 * a)
  if (x1 > 0): 
    return x1
  else:
    return (-b - sqrt_val)/(2 * a)
for i in range(0, n):
  t1 = t0 + h
  J1 = quadric(-11/90, -1.1 + 1/9*(R0 + J0), J0)
  R1 = (R0 + J0 -1.1*J1) / 0.9
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

y0 = [2, 4] #vector [R, J]
n = 10
t = np.linspace(0, n, int(n/0.1))
def pend(y, t,):
    R, J = y
    dydt = [ R*(1 - J), J*(R - 1)]
    return dydt



#plot solution
sol = odeint(pend, y0, t)
plt.plot(t, sol[:, 0], 'b', label='R(t)')
plt.plot(t, sol[:, 1], 'g', label='J(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.title("Exact Solution (with help of scipy library in Python")
plt.show()