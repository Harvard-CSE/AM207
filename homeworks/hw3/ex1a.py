#!/usr/bin/env python3
import korali
import math
import numpy as np 

def negative_ackley(p):
    ''' this function returns -f where f is the function given in the homework sheet 
        korali will maximise -f i.e minimise f. 
    '''
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
    r1 = a*np.exp(-b*np.sqrt(sum1))
    r2 = np.exp(sum2)

    p["F(x)"] = r1 + r2 - a - np.exp(1)
    
    grad = [0.]*dim
    for i in range(dim):
      grad[i] = r1*-1*b*0.5/np.sqrt(sum1)*1.0/dim*2.0*x[i]
      grad[i] -= r2*1.0/dim*np.sin(c*x[i])*c

    p["Gradient"] = grad

k = korali.Engine()
e = korali.Experiment()

# Configuring Problem
e["Random Seed"] = 0xC0FEE 
e["Problem"]["Objective Function"] = negative_ackley
e["Problem"]["Type"] = "Optimization"

dim = 2

# Defining the problem's variables.
for i in range(dim):
    e["Variables"][i]["Name"] = "X" + str(i)
    e["Variables"][i]["Initial Value"] = 5
    e["Variables"][i]["Lower Bound"] = -32.768
    e["Variables"][i]["Upper Bound"] = 32.768
    e["Variables"][i]["Initial Standard Deviation"] = 5
    e["Variables"][i]["Initial Mean"] = 5
    

# Configuring CMA-ES parameters
e["Solver"]["Type"] = "Optimizer/CMAES"
e["Solver"]["Population Size"] = 32
e["Solver"]["Mu Value"] = 8

e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-32
e["Solver"]["Termination Criteria"]["Max Generations"] = 200
    
# Configuring results path
e["File Output"]["Enabled"] = True
e["File Output"]["Path"] = '_korali_result_cmaes'
e["File Output"]["Frequency"] = 1

# Running Korali
k.run(e)
    