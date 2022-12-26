#!/usr/bin/env python
import numpy as np


# multi dimensional problem (ackley)
def negative_ackley(p):
    x = p["Parameters"]
    a = 20.
    b = 0.2
    c = 2.*np.pi
    dim = len(x)

    sum1 = 0.
    sum2 = 0.
    for i in range(dim):
        sum1 += x[i]*x[i]
        sum2 += np.cos(c*x[i])

    sum1 /= dim
    sum2 /= dim
    r1 = -a*np.exp(-b*np.sqrt(sum1))
    r2 = -np.exp(sum2)

    p["F(x)"] = r1 + r2 + a + np.exp(1)
    
    grad = [0.]*dim
    for i in range(dim):
      grad[i] = r1*-1*b*0.5/np.sqrt(sum1)*1.0/dim*2.0*x[i]
      grad[i] -= r2*1.0/dim*np.sin(c*x[i])*c

    p["Gradient"] = grad