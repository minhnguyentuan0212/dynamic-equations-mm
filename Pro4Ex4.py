import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import math
#IMPLICIT METHOD SOLUTIONt
#configuration
R = [1]
J = [1]
t = [0]
n = 100
#initial value
R0 = 1
J0 = 1
h = 0.1
t0 = 0

#NewtonRaphson
def NewtonRaphson(initR, initJ, R, J, t, delta = 0.01):
  initDelta = np.array([[initR], [initJ]])
  RJ = np.array([[R], [J]])
  deltaRJ = np.array([[1000.], [1000.]])
  while abs(deltaRJ.item(0)) > delta or abs(deltaRJ.item(1)) > delta:
    H = np.array([[0.01*math.e**RJ.item(0) + 0.1*RJ.item(1) + initDelta.item(0) - RJ.item(0)], [-0.005*math.e**(2*RJ.item(1)) - 0.3*RJ.item(0) + initDelta.item(1) - RJ.item(1)]])
    J = np.array([[0.01*math.e**RJ.item(0) - 1, 0.1], [-0.3, -0.01*math.e**(2*RJ.item(1)) - 1]])
    deltaRJ = np.matmul(np.linalg.inv(J), H)
    RJ = np.subtract(RJ, deltaRJ)
  return RJ
#iteration
for i in range(0, n):
  t1 = t0 + h
  #processing
  RJ = NewtonRaphson(R0, J0, R0, J0, t1)
  R1 = RJ.item(0)
  J1 = RJ.item(1)
  #assign 
  R.append(R1)
  J.append(J1)
  t.append(t0)
  R0 = R1
  J0 = J1
  t0 = t1
#plot solution
plt.plot(t, R, 'b', label = "R line")
plt.plot(t, J, 'g', label = "J line")
plt.legend()
plt.legend(loc='best')
plt.grid()
plt.title("Implicit Euler solution")
plt.show()
#EXACT SOLUTION 

y0 = [1, 1] #vector [R, J]
n = 10
t = np.linspace(0, n, int(n/0.1))
def pend(y, t,):
    R, J = y
    dydt = [ 0.1*math.e**(R) + J ,  -0.05*math.e**(2*J) - 3*R]
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