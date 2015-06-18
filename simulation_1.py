#-*- encoding:utf-8 -*-
__author__ = 'qsong'

from scipy.stats import poisson
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import math
import logging
import os
import random

import matplotlib.legend_handler as HandlerLine2D

# 以下参数定义为全局变量
# all numerical results take the following values
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

# CDMA related function
def average_power(intensity, epsilon, mu, g, l, tau_s, w, p_f, r_0, r_i):
    Nc = cal_code_length(intensity, epsilon, p_f)
    n_mean = cal_n_bar(intensity, epsilon)
    mu_t = pow(2, l*Nc/w/tau_s)-1
    denominator = mu*g*(Nc/mu_t-n_mean+1)
    avg_power_dbm = 10*math.log(1000*1.0/denominator, 10)
    print 'mu_t', mu_t
    print intensity, '-->', n_mean, '-->', Nc, '-->', 1.0/denominator, '-->', avg_power_dbm
    return avg_power_dbm

def average_power2(intensity, epsilon, mu, g, l, tau_s, w, p_f, r_0, r_i):
    Nc = cal_code_length(intensity, epsilon, p_f)
    n_mean = cal_n_bar(intensity, epsilon)
    mu_t = pow(2, l*Nc/w/tau_s)-1
    denominator = mu*(Nc/mu_t-n_mean+1)
    avg_power_dbm = 10*math.log(1000*15.0/denominator, 10)
    print 'mu_t', mu_t
    print intensity, '-->', n_mean, '-->', Nc, '-->', 1.0/denominator, '-->', avg_power_dbm
    return avg_power_dbm

def average_energy():
    pass
    # return  average_power()*tau_s/l

def cal_n_bar(intensity, epsilon):
    """
        This function is used to calculate the (1-epsilon)th percentile of arrivals
    """
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    k = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[k]
        k += 1
    return k


def cal_code_length(intensity, epsilon, p_f):
    """
        This function is used to calculate CDMA code length
    """
    p_c = (p_f-epsilon)/(1.0-epsilon)
    n_mean = cal_n_bar(intensity, epsilon)
    return math.ceil(math.log(1+(n_mean-1)/p_c, 2))


def CDMA_average_with_variable_SNR(intensity, epsilon, p_f, i=1):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, i*L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    p_t = CONSTANT/(MU*(N_c/mu_t-n_mean+1))
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}, N_c is :{2}, N_bar:{3}".format(p_t, i, N_c, n_mean)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    return p_t

def CDMA_average_with_mean_SNR(intensity, epsilon, p_f, prob_dist):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, i*L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    # p_t = CONSTANT/(MU*(N_c/mu_t-n_mean+1))
    p_t = CONSTANT*sum([q[1]/(MU*(N_c/(pow(2, q[0]*L*N_c/W/TAU_S)-1)-n_mean+1)) for q in prob_dist])
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}, N_c is :{2}, N_bar:{3}".format(p_t, i, N_c, n_mean)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    return p_t

def CDMA_ratio_target_SNR(intensity, epsilon, p_f, i=1):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, i*L*N_c/W/TAU_S)-1
    # The maximum target SNR can be N_c/(n_mean-1)
    # The minimum target SNR is pow(2, i*L*N_c/W/TAU_S)-1, namely mu_t
    return (N_c/(n_mean-1))/mu_t

def CDMA_average_power_variant_g(intensity, epsilon, p_f, bound, i=1):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, i*L*N_c/W/TAU_S)-1
    CONSTANT_1 = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    CONSTANT_2 = 1.0*(math.log(1+bound[1], np.e)-math.log(1+bound[0], np.e))/(bound[1]-bound[0])
    p_t = CONSTANT_2*CONSTANT_1/(MU*(N_c/mu_t-n_mean+1))
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}".format(p_t, i)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    return p_t


# last time to modify this function:
# The transmit power is function of three independent random variables, distance r, power control error x and a
# log normal y.

def CDMA_average_power_normal_g(intensity, epsilon, p_f, bound, i=1, mu=0, delta=1):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, i*L*N_c/W/TAU_S)-1
    CONSTANT_1 = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

    # def h(x, mu, delta):
    #     return 1.0/(math.sqrt(2*math.pi*delta**2))*math.exp(-(x-mu)**2/(2*delta**2))*(10**(x/10.0))
    #
    def h(x, mu, delta):
        return 1.0/(math.sqrt(2*math.pi*delta**2))*math.exp(-(x-mu)**2/(2*delta**2))/(N_c*mu_t*10**(x/10.0)-(n_mean-1))

    # def h(y, x, mu, sigma):
    #     E = 0.2303
    #     var_Y = math.log(1+(np.e**((E*sigma)**2)-1)/(n_mean-1), np.e)
    #     sigma_Y = math.sqrt(var_Y)
    #     mu_Y = math.log(n_mean-1, np.e) + 0.5*((E*sigma)**2 - var_Y)
    #     denom = 2*math.pi*sigma*sigma_Y*y*((N_c/mu_t)*10**(x/10.0)-y)
    #     # print "y", y
    #     nom = np.e**(-x**2/(2*sigma**2)-(math.log(y, np.e)-mu_Y)**2/(2*var_Y))
    #     result = 0
    #     try:
    #         result = nom/denom
    #     except ZeroDivisionError:
    #         print "x,y:", x, y
    #         result = nom/(denom+0.001)
    #     return result


    CONSTANT_2 = integrate.quad(h, bound[0], bound[1], args=(mu, delta))[0]

    # CONSTANT_2 = integrate.dblquad(
    #     h, -100, 100, lambda x:0.01, lambda x:100, args=(mu, delta)
    # )[0]

    p_t = CONSTANT_2*CONSTANT_1/MU

    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1} with intensity {2}".format(p_t, i, intensity)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    return p_t


def MT_CDMA_average_power_normal_g(intensity, epsilon, p_f, bound, i=1, mu=0, delta=1):

    p_c = (p_f-epsilon)/(1.0-epsilon)

    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]

    cumulative = 0.0

    n_mean = 0

    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1

    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))

    mu_t = pow(2, i*L*N_c/W/TAU_S)-1

    CONSTANT_1 = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

    CONSTANT_2 = [
        N_c/mu_t*10**(random.gauss(mu, delta)/10) - sum([10**(random.gauss(mu, delta)/10.0) for i in range(int(n_mean-1))])
        for i in range(10000)
    ]

    # CONSTANT_2 = np.mean(CONSTANT_2)

    D = [CONSTANT_1/element/MU for element in CONSTANT_2]



    p_t = np.mean(D)

    try:

        p_t = 10*math.log(1000*p_t, 10)

    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1} with intensity {2}".format(p_t, i, intensity)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))

    return p_t

def average_power4(intensity, epsilon, p_f):
    p_c = (p_f-epsilon)/(1.0-epsilon)
    n_mean = intensity
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    p_t = CONSTANT/(MU*(N_c/mu_t-n_mean+1))
    p_t = 10*math.log(1000*p_t, 10)
    logger.debug("lamdba:{0:5s} N_bar:{1:4s} N_C:{2:4s} p_t:{3:4s}".format(intensity, n_mean, N_c, p_t))
    return p_t



# FDMA related functions

def FDMA_average_power(intensity, epsilon):
    """
        This function is used to reproduce the numerical results given in Dhillion's Paper...
    """
    # p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    K = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[K]
        K += 1
    # N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    # mu_t = pow(2, L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    p_t = (2**(K*L/W/TAU_S)-1.0)/MU/K*CONSTANT
    p_t = 10*math.log(1000*p_t, 10)
    logger.debug("lamdba:{0} K:{1} p_t:{2}".format(intensity, K, p_t))
    return p_t


def VL_FDMA_average_power(intensity, epsilon):
    """
        This function is used to explore the impact when packet length is a variable instead of a constant
    """
    # p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    K = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[K]
        K += 1
    # N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    # mu_t = pow(2, L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    p_t = (2**(5*K*L/W/TAU_S)-1.0)/MU/K*CONSTANT
    p_t = 10*math.log(1000*p_t, 10)
    logger.debug("lamdba:{0} K:{1} p_t:{2}".format(intensity, K, p_t))
    return p_t


# SIC related functions
# 虽然我一直不太明白，这个所谓的SIC倒是是个什么鸟玩意儿。。。
def SIC_average_power(intensity, epsilon):
    # p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    K = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[K]
        K += 1
    # N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    # mu_t = pow(2, L*N_c/W/TAU_S)-1
    CONSTANT = (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    power_list = [2**((k-1)*L/W/TAU_S)*(2**(L/W*TAU_S)-1)/MU*CONSTANT for k in range(1, K+1)]
    p_t = sum(power_list)/len(power_list)
    # p_t = 2**((K-1)*L/W/TAU_S)*(2**(L/W*TAU_S)-1)/MU*CONSTANT
    p_t = 10*math.log(1000*p_t, 10)
    logger.debug("lamdba:{0} K:{1}, p_t:{2}".format(intensity, K, p_t))
    return p_t


# 添加新方法，用以评估引入power control error之后，CDMA条件下，average transmit power
def lower_bound_imperct_cdmd(intensity, epsilon, p_f, mu=0, delta=1):
    """
        经过数学分析，in case of imperfect CDMA, the upper bound of average transmit power is:
            \mu_t \exp{1/2 (1n(1/10)\sigma/10)^2} E[1/(mu*g)]
        其中，经过计算：
        E[1/(mu g)] = mu * (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

        @:parameter
        indensity: float, arrival request density
        mu: float, default value 0, mean of normal distribution
        delta: float, default value 1, variance of normal distribution
    """
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, L*N_c/W/TAU_S)-1
    print "mu_t", mu_t, "N_c", N_c, "N", n_mean, "Nc - mu_t*(N-1)", N_c-mu_t*(n_mean-1)

    CONSTANT_1 = (1.0/N_c)*(1.0/MU)*(RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    CONSTANT_2 = math.exp(0.5*(mu+math.log(0.1, np.e)*delta/10.0)**2)
    p_t = mu_t * CONSTANT_1*CONSTANT_2
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}".format(p_t, i)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    print "CONSTANT_1,CONSTANT_2,p_t", CONSTANT_1, CONSTANT_2, p_t
    return p_t


# 添加新方法，用以评估引入power control error之后，CDMA条件下，average transmit power
def upper_bound_imperct_cdmd(intensity, epsilon, p_f, mu=0, delta=1):
    """
        经过数学分析，in case of imperfect CDMA, the upper bound of average transmit power is:
            \mu_t \exp{1/2 (1n(1/10)\sigma/10)^2} E[1/(mu*g)]
        其中，经过计算：
        E[1/(mu g)] = mu * (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

        @:parameter
        indensity: float, arrival request density
        mu: float, default value 0, mean of normal distribution
        delta: float, default value 1, variance of normal distribution
    """
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, L*N_c/W/TAU_S)-1

    CONSTANT_1 = (1.0/MU)*(RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    CONSTANT_2 = (n_mean-1)*math.exp(0.5*(mu+math.log(10, np.e)*delta/10.0)**2) + CONSTANT_1
    CONSTANT_3 = math.exp(0.5*(mu+math.log(0.1, np.e)*delta/10.0)**2)
    p_t = (mu_t/N_c) * CONSTANT_3 * CONSTANT_2
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}".format(p_t, i)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    print "CONSTANT_1,CONSTANT_2, CONSTANT_3=", CONSTANT_1, CONSTANT_2, CONSTANT_3
    # 如果p_t的值小于1，则返回p_t.不然返回1
    return p_t
        # if p_t < 30.0 else 30.0


# 添加新方法，用以评估引入power control error之后，CDMA条件下，average transmit power
def mean_imperct_cdmd(intensity, epsilon, p_f, mu=0, delta=1):
    """
        经过数学分析，in case of imperfect CDMA, the upper bound of average transmit power is:
            \mu_t \exp{1/2 (1n(1/10)\sigma/10)^2} E[1/(mu*g)]
        其中，经过计算：
        E[1/(mu g)] = mu * (RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)

        @:parameter
        indensity: float, arrival request density
        mu: float, default value 0, mean of normal distribution
        delta: float, default value 1, variance of normal distribution
    """
    p_c = (p_f-epsilon)/(1.0-epsilon)
    proba_distribution_poisson = [poisson.pmf(k, intensity) for k in range(int(intensity+200))]
    cumulative = 0.0
    n_mean = 0
    while cumulative <= 1-epsilon:
        cumulative += proba_distribution_poisson[n_mean]
        n_mean += 1
    N_c = math.ceil(math.log(1+(n_mean-1)/p_c, 2))
    mu_t = pow(2, L*N_c/W/TAU_S)-1

    CONSTANT_1 = (1.0/MU)*(RADIUS**(GAMMA+2)-RADIUS_INNER**(GAMMA+2))*2/((RADIUS**2-RADIUS_INNER**2)*(GAMMA+2)*RADIUS**GAMMA)
    CONSTANT_2 = (n_mean-1)*math.exp(0.5*(mu+math.log(10, np.e)*delta/10.0)**2)*(N_c-mu_t*(n_mean-1))
    # CONSTANT_3 = math.exp(0.5*(mu+math.log(0.1, np.e)*delta/10.0)**2)
    p_t = (mu_t) * CONSTANT_1 / CONSTANT_2
    try:
        p_t = 10*math.log(1000*p_t, 10)
    except ValueError:
        print "Math domain error,p_t is:{0} when i is:{1}".format(p_t, i)
        print "CONSTANT 2:{0}".format(CONSTANT_2)

    logger.debug("lamdba:{0} N_bar:{1} N_C:{2} p_t:{3}".format(intensity, n_mean, N_c, p_t))
    print "CONSTANT_1,CONSTANT_2", CONSTANT_1, CONSTANT_2
    # 如果p_t的值小于1，则返回p_t.不然返回1
    return p_t
        # if p_t < 30.0 else 30.0



if __name__ == '__main__':


    EPSILON = 0.01
    P_F = 0.06
    DELTA = 0.0
    P_C = (P_F-EPSILON)/(1-EPSILON)

    MAX_LAMBDA = 420

    FIGURES_FOLDED = '/Users/qsong/Documents/PhD_Doc/publication/Figures'

    # 利用Python logging模块，记录脚本运行信息，有利于debug...
    logging.basicConfig(
        filename=os.path.join(os.getcwd(), 'debug_log.txt'),
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s - %(levelname)s: %(message)s'
    )
    logger = logging.getLogger(__name__)
    logger.debug("General system parameter:")
    logger.debug("The maximum transmit power: {0}".format(P_MAX))
    logger.debug("The outer radii: {0}".format(RADIUS))
    logger.debug("The inner radii: {0}".format(RADIUS_INNER))
    logger.debug("The path loss exponent: {0}".format(GAMMA))
    logger.debug("The reference SNR: {0}".format(MU))
    logger.debug("The total bandwidth: {0}".format(W))
    logger.debug("The length of unit time slot: {0}\n".format(TAU_S))
    logger.debug("CDMA-related system parameter:")
    logger.debug("The overall failure probability in CDMA random access: {0}".format(P_F))
    logger.debug("The outage probability: {0}".format(DELTA))
    logger.debug("The 1-epsilon th percentile of arrivals: {0}".format(EPSILON))

    # figure 1, plot the reference figure, presented in Dhillon's paper
    plt.figure(1, figsize=(12, 4.5), dpi=120)
    x_BS_load = np.linspace(1.1, 1300, num=100)
    y_avg_power_1_1 = [CDMA_average_with_variable_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F) for intensity in x_BS_load]
    line_cdma_avg_power, = plt.plot(x_BS_load, y_avg_power_1_1, linestyle='--', marker='s', markevery=4, color='b', label='CDMA:Mean')
    logger.debug("-----------------------------------------FDMA---------------------------------------------")
    y_avg_power_2 = [FDMA_average_power(intensity=intensity, epsilon=EPSILON) for intensity in x_BS_load]
    logger.debug("-----------------------------------------SIC---------------------------------------------")
    y_avg_power_3 = [SIC_average_power(intensity=intensity, epsilon=EPSILON) for intensity in x_BS_load]
    line_fdma_avg_power, =plt.plot(x_BS_load, y_avg_power_2, linestyle='--', marker='D', markevery=4, color='r', label='FDMA: Mean')
    line_sic_avg_power, = plt.plot(x_BS_load, y_avg_power_3, linestyle='--', marker='v', markevery=4, color='m', label='SIC: Mean')
    plt.xlabel(r'BS load ($\lambda$ arrivals per second)', fontsize=12)
    plt.ylabel(r'Average transmit power (in dBm)', fontsize=12)
    plt.xlim([0, 1300])
    # plt.ylim([-4, 16])
    plt.grid()
    first_legend = plt.legend(
        handles=[
            line_cdma_avg_power,
            line_fdma_avg_power,
            line_sic_avg_power
        ],

        loc=2
    )
    ax = plt.gca().add_artist(first_legend)
    plt.savefig(os.path.join(FIGURES_FOLDED, "Reproduction_Dhillon.png"))


    # figure 2, plot the relationship between BS load and target
    # ------------------------------------------------------------------------------------------------------------------
    plt.figure(2, figsize=(12, 4.5), dpi=120)
    x_BS_load_2 = np.linspace(100, 1500, num=100)
    y_tgt_SNR_ratio = [CDMA_ratio_target_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F) for intensity in x_BS_load_2]
    plt.xlabel(r'BS load ($\lambda$ arrivals per second)', fontsize=12)
    plt.ylabel(r'Ratio of $\mu_{t_{max}}/\mu_{t_{min}}$', fontsize=12)
    plt.plot(x_BS_load_2, y_tgt_SNR_ratio, linestyle='-', marker='*', markevery=4, color='g')
    plt.grid()
    plt.xlim([100, 1500])
    plt.yticks(range(0, 13), [str(i) for i in range(0, 13)])
    plt.xticks(range(100, 1600, 100), [str(i) for i in range(100, 1600, 100)])
    plt.annotate('(1400, 1)', xy=(1400, 1.1), xytext=(1400, 3.1), arrowprops=dict(arrowstyle="->"))
    plt.annotate('(800, 1.64)', xy=(800, 1.64), xytext=(800, 3.1), arrowprops=dict(arrowstyle="->"))
    3.06
    plt.annotate('(420, 3.06)', xy=(420, 3.1), xytext=(420, 4.1), arrowprops=dict(arrowstyle="->"))
    plt.savefig(os.path.join(FIGURES_FOLDED, "Ratio-target-SNR-max-min.png"))



    # figure 3, plot the evalution when playing with target SNR mu_t (packet length)
    # ------------------------------------------------------------------------------------------------------------------
    plt.figure(3, figsize=(12, 4.5), dpi=120)
    x_BS_load = np.linspace(1.1, MAX_LAMBDA, num=100)
    y_avg_power_1_1 = [CDMA_average_with_variable_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F) for intensity in x_BS_load]
    # 极限情况： cell中所有的packet均为3L，此时为 Average transmit power曲线之上限
    y_avg_power_1_3 = [CDMA_average_with_variable_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F, i=3) for intensity in x_BS_load]
    y_avg_power_1_5 = [CDMA_average_with_mean_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F, prob_dist=[(0.8, 0.2), (1.0, 0.15), (1.2, 0.15), (1.6, 0.1), (2.0, 0.1), (2.5, 0.15), (3.0, 0.15)]) for intensity in x_BS_load]
    # FDMA 的参考线
    y_avg_power_2_1 = [FDMA_average_power(intensity=intensity, epsilon=EPSILON) for intensity in x_BS_load]
    y_avg_power_2_2 = [VL_FDMA_average_power(intensity=intensity, epsilon=EPSILON) for intensity in x_BS_load]



    line_cdma_avg_power, = plt.plot(x_BS_load, y_avg_power_1_1, linestyle='--', marker='s', markevery=4, color='b', label='CDMA: reference case(constant packet length)')
    line_cdma_avg_power_variant_mut, = plt.plot(x_BS_load, y_avg_power_1_3, linestyle='--', marker='+', markevery=4, color='k', label='CDMA: limit case where constant packet are 3L')
    line_cdma_avg_power_mean_mut, = plt.plot(x_BS_load, y_avg_power_1_5, linestyle='--', marker='*', markevery=4, color='g', label='CDMA: general case where packet length is discret random variable')


    line_fdma_avg_power, =plt.plot(x_BS_load, y_avg_power_2_1, linestyle='--', marker='D', markevery=4, color='r', label='FDMA: Mean')
    line_fdma_avg_power_1, =plt.plot(x_BS_load, y_avg_power_2_2, linestyle='--', marker='o', markevery=4, color='m', label='FDMA: with 3L')

    plt.grid()
    plt.xlim([0, MAX_LAMBDA])
    # plt.legend(handler_map={line_cdma_avg_power: HandlerLine2D(numpoints=4)})
    third_legend = plt.legend(
        handles=[
            line_cdma_avg_power,
            line_cdma_avg_power_variant_mut,
            line_cdma_avg_power_mean_mut,
            line_fdma_avg_power,
            line_fdma_avg_power_1
        ],
        loc=2
    )
    ax = plt.gca().add_artist(third_legend)
    plt.xlabel(r'BS load ($\lambda$ arrivals per second)', fontsize=12)
    plt.ylabel(r'Average transmit power (in dBm)', fontsize=12)
    plt.savefig(os.path.join(FIGURES_FOLDED, "Average-Transmit-Power-Variant-Target-SNR.png"))


    # # change: 2015-06-16 注释掉关于 channel g
    # # figure 4, plot the evaluation when playing with channel gain g
    # # ------------------------------------------------------------------------------------------------------------------
    # plt.figure(4, figsize=(12, 4.5), dpi=120)
    # x_BS_load = np.linspace(1.1, MAX_LAMBDA, num=100)
    # y_avg_power_1_1 = [CDMA_average_with_variable_SNR(intensity=intensity, epsilon=EPSILON, p_f=P_F) for intensity in x_BS_load]
    # line_cdma_avg_power, = plt.plot(x_BS_load, y_avg_power_1_1, linestyle='--', marker='s', markevery=4, color='b', label='CDMA:reference')
    # #y_avg_power_variant_g_1 = [CDMA_average_power_variant_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-0.75, 0.75)) for intensity in x_BS_load]
    # #y_avg_power_variant_g_2 = [CDMA_average_power_variant_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-0.25, 0.25)) for intensity in x_BS_load]
    # #y_avg_power_variant_g_3 = [CDMA_average_power_variant_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-0.5, 0.5)) for intensity in x_BS_load]
    # y_avg_power_variant_g_1 = [MT_CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=1) for intensity in x_BS_load]
    #
    # # y_avg_power_variant_g_2 = [CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=2) for intensity in x_BS_load]
    # #
    # # y_avg_power_variant_g_3 = [CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=3) for intensity in x_BS_load]
    # # y_avg_power_variant_g_4 = [CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=4) for intensity in x_BS_load]
    # #
    # # y_avg_power_variant_g_4 = [CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=4) for intensity in x_BS_load]
    # for i in range(2):
    #     plt.scatter(
    #         x_BS_load,
    #         [MT_CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=1) for intensity in x_BS_load],
    #         color='red'
    #     )
    #
    # for i in range(2):
    #     plt.scatter(
    #         x_BS_load,
    #         [MT_CDMA_average_power_normal_g(intensity=intensity, epsilon=EPSILON, p_f=P_F, bound=(-100, 1000), delta=2) for intensity in x_BS_load],
    #         color='green'
    #     )
    # # line_cdma_avg_power_variant_g_2, = plt.plot(x_BS_load, y_avg_power_variant_g_2, linestyle='--', marker='*', markevery=4, color='y', label='CDMA:$\sigma=2$')
    # # line_cdma_avg_power_variant_g_3, = plt.plot(x_BS_load, y_avg_power_variant_g_3, linestyle='--', marker='*', markevery=4, color='c', label='CDMA:$\sigma=3$')
    # # line_cdma_avg_power_variant_g_4, = plt.plot(x_BS_load, y_avg_power_variant_g_4, linestyle='-', marker='*', markevery=4, color='c', label='CDMA:$\sigma=4$')
    #
    # plt.grid()
    # plt.xlim([0, MAX_LAMBDA])
    # # # plt.legend(handler_map={line_cdma_avg_power: HandlerLine2D(numpoints=4)})
    # # first_legend = plt.legend(
    # #     handles=[
    # #         line_cdma_avg_power,
    # #         line_cdma_avg_power_variant_g_1,
    # #         line_cdma_avg_power_variant_g_2,
    # #         line_cdma_avg_power_variant_g_3,
    # #         line_cdma_avg_power_variant_g_4
    # #     ],
    # #     loc=2
    # # )
    # # ax = plt.gca().add_artist(first_legend)
    # plt.xlabel(r'BS load ($\lambda$ arrivals per second)', fontsize=12)
    # plt.ylabel(r'Average transmit power (in dBm)', fontsize=12)
    # plt.savefig(os.path.join(FIGURES_FOLDED, "Average-Transmit-Power-Variant-Channel-Gain.png"))

    plt.show()