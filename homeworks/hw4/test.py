import korali
import numpy as np
import pandas as pd


df = pd.read_csv('data.csv')
x = df["x"].to_numpy()
# y = df["y"].to_numpy()
n = len(x)

def log_likelihood(ks):
    a, b = ks["Parameters"]
    log_like = np.sum(np.log(p_x(a,b)))
    ks["logLikelihood"] =  log_like
    

def my_log_likelihood(ks):
    a,b =  ks["Parameters"]
    log_like = np.sum(np.log(p_x(a,b)))
    ks["logLikelihood"] = log_like

def p_x(a,b):
    'pointwise likelihood'
    return (1/np.pi) * ((b)/(b**2 + (a - x)**2))


e = korali.Experiment()

e["Problem"]["Type"] = "Bayesian/Custom"
e["Problem"]["Likelihood Model"] = my_log_likelihood

e["Solver"]["Type"] = "Sampler/TMCMC"
e["Solver"]["Population Size"] = 5000
e["Solver"]["Target Coefficient Of Variation"] = 1.0
e["Solver"]["Covariance Scaling"] = 0.2

e["Distributions"] = [
    {"Name": "Prior a",
        "Type": "Univariate/Uniform",
        "Minimum": np.max(x) - np.std(x),
        "Maximum": np.max(x) + np.std(x)},

    {"Name": "Prior b",
        "Type": "Univariate/Uniform",
        "Minimum": np.max(x) - np.mean(x) - 2*np.std(x),
        "Maximum": np.max(x) - np.mean(x) + 2*np.std(x)}
]


e["Variables"] = [
    {"Name": "a", "Prior Distribution": "Prior a"},
    {"Name": "b", "Prior Distribution": "Prior b"}
]

e["Store Sample Information"] = True

e["Console Output"]["Verbosity"] = "Detailed"
k = korali.Engine()
k.run(e)