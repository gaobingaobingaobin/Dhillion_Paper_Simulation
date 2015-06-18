import random
import math
import numpy as np

N_c = 13
mu_t = 2**(N_c/1000.0) - 1
n_mean = 353


P_MAX = 1.0
RADIUS = 1000.0
# In Dhillon's paper, no give the value of inner radius. We suppose it as 200 meters.
RADIUS_INNER = 50.0
GAMMA = 3.0
# Be careful, the unit of MU is decibel.
# In decibel manner, SNR_dB = 10*log10(SNR) = > SNR = 10^(SNR_db/10)
MU = pow(10, -3.0/10)
W = 1e6
TAU_S = 1.0
L = 1e3

rand1 = random.gauss(0, 1)/10
p_error1 = 10**(rand1)
intefer = sum([10**(random.gauss(0,1)/10.0) for i in range(n_mean-1)])

CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

mu = 0
sigma = 2
C = [N_c/mu_t*10**(random.gauss(mu, sigma)/10) - sum([10**(random.gauss(mu, sigma)/10.0) for i in range(n_mean-1)])  for i in range(10000)]
D = [CONSTANT/element/MU for element in C]

P_t = 10*math.log(1000*np.mean(D), 10)

print P_t