import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import math

sys.path.append("../src/model_A")

import runnerA
import fitness_functions_model_A as ffmA

sys.path.append("../src/model_B")

import runnerB
import fitness_functions_model_B as ffmB

sys.path.append("../src/model_C")

import runnerC
import fitness_functions_model_C as ffmC

solA = runnerA.PySolver()

solB = runnerB.PySolver()

solC = runnerC.PySolver()

b = 2.0
q = 0.1
d = 0.1
alpha = 0.5
gamma = 0.5

delta = 0.1
eta = 0.9

R0_upper_limit = 10
labels = [["(A)", "(B)", "(C)"], ["(D)", "(E)", "(F)"], ["(G)", "(H)", "(I)"]]

c1 = 1.0

c2 = -3.0
seed = 31

resolution = 300


S_density = I_F_density = 3.0
I_Q_density = 0.0
trait_space_length = 101

evo_steps = 2000

rho_max = 1.0
beta_max = 1.0

depth_max = 3

scenario = 3

rhos = [i/resolution for i in range(resolution + 1)]
sigma_extinctions = []
for rho in rhos:
    beta_val = ffmA.beta_func(beta_max, rho, c1, c2)
    if rho != 1:
        extinction_sigma = ((b - d) * beta_val - b * q * (d + alpha + gamma)) / (beta_val * eta *(delta - 1) *(-1 + rho)* (b - d))
    else:
        if (b - d) * beta_val - b * q * (d + alpha + gamma) > 0:
            extinction_sigma = 1
        else:
            extinction_sigma = 0
    sigma_extinctions.append(extinction_sigma)
sigma_extinction = max(sigma_extinctions)
print(sigma_extinction)


sigma_calc = ffmA.dbetadrho_func(beta_max, 0, c1, c2)/(eta*(delta - 1)*(ffmA.beta_func(beta_max, 0, c1, c2) - ffmA.dbetadrho_func(beta_max, 0, c1, c2)))
print(sigma_calc)
evolutionary_stabilities = []
convergence_stabilities = []
evolved_rhos = []
evolved_sigmas = []
if sigma_extinction < 1:
    
    sigmas = [sigma_extinction*i/resolution for i in range(resolution + 1)]
#     print(sigmas)
else:
    sigmas = [1.0*i/resolution for i in range(resolution + 1)]
    
# flop
for sigma in sigmas:
    attractors, repellers = ffmA.find_sing_strats(solA, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, sigma, seed, scenario)
    for rho in attractors:
        evolved_rhos.append(rho)
        evolved_sigmas.append(sigma)
        
# print(evolved_rhos)
fig = plt.figure(figsize = (15, 10))
gs = fig.add_gridspec(3,3)
ax = [[fig.add_subplot(gs[j,i]) for i in range(3)] for j in range(3)]

prevs, prevs_no_evos = ffmA.calculate_disease_prevalance(solA, evolved_rhos, evolved_sigmas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)
R0s, R0s_no_evo = ffmA.R0s_model_A(solA, evolved_rhos, evolved_sigmas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)


ax[0][0].plot(evolved_sigmas, evolved_rhos)
# plt.scatter(evolved_sigmas, simulated_rhos)
ax[0][0].plot([sigma_calc, sigma_calc], [0, 1.1], "k--")
if sigma_extinction < 1:
    ax[0][0].plot([sigma_extinction, sigma_extinction], [0, 1.1], "k-.")
ax[0][0].set_ylim([0,1.1])
ax[0][0].set_xlim([0,1.0])
ax[0][0].set_ylabel(r"Evolved false negativity, $\rho^*$")
ax[0][0].set_title("Model A")
ax[0][0].text(0.03, 0.90,labels[0][0], transform=ax[0][0].transAxes)

ax[1][0].plot(evolved_sigmas, prevs)
ax[1][0].plot(evolved_sigmas, prevs_no_evos)
ax[1][0].plot([sigma_calc, sigma_calc], [0, 1.1], "k--")
if sigma_extinction < 1:
    ax[1][0].plot([sigma_extinction, sigma_extinction], [0, 1.1], "k-.")
ax[1][0].set_ylim([0,0.2])
ax[1][0].set_xlim([0,1])
ax[1][0].set_ylabel("Pathogen Prevalance")
# ax[1][0].set_xlabel(r"Testing probability, $\sigma$")
ax[1][0].text(0.03, 0.90,labels[1][0], transform=ax[1][0].transAxes)

ax[2][0].plot(evolved_sigmas, R0s)
ax[2][0].plot(evolved_sigmas, R0s_no_evo)
ax[2][0].plot([sigma_calc, sigma_calc], [0, R0_upper_limit], "k--")
# ax[2].plot([sigma_calc, sigma_calc], [0, 1.1], "k--")
ax[2][0].plot([0,1],[1,1], "k--")
if sigma_extinction < 1:
    ax[2][0].plot([sigma_extinction, sigma_extinction], [0, R0_upper_limit], "k-.")
# ax[2].set_ylim([0,1.1])
ax[2][0].set_xlim([0,1])
ax[2][0].set_ylim([0,R0_upper_limit])
ax[2][0].set_ylabel(r"Basic Reproduction Number, $R_0$")
ax[2][0].set_xlabel(r"Testing probability, $\sigma$")
ax[2][0].text(0.03, 0.90,labels[2][0], transform=ax[2][0].transAxes)
# plt.savefig("../outputs/test_sigma_variation.png", bbox_inches = "tight")
# plt.close()

d_tilde = d + alpha + gamma
if delta == 0 and eta == 1.0:
    zeta_calc = (-ffmB.dbetadrho_func(beta_max,0,c1,c2)*(d + alpha + gamma))/(ffmB.beta_func(beta_max,0,c1,c2) + ffmB.dbetadrho_func(beta_max,0,c1,c2))
    
    zeta_limits = []
    zeta_res = 300
    for i in range(zeta_res + 1):
        rho = rho_max*i/zeta_res
        beta = ffmB.beta_func(beta_max,rho,c1,c2)
#         print(beta)
        if rho != 1:
#             temp_zeta = ((beta*(b - d)/(b*q)) - (d + alpha + gamma))/(1 - rho)
            temp_zeta = (((d + alpha + gamma)*q - beta)*b + beta*d)/(b*q*(rho - 1))
        else:
            temp_zeta = beta*(b - d)/(b*q*(d + alpha + gamma))
        zeta_limits.append(temp_zeta)
    zeta_limit = max(zeta_limits)
#     print(zeta_limit, zeta_limits)
    zetas = [zeta_limit*i/resolution for i in range(resolution + 1)]
else:
    zeta_divider = 2*ffmB.dbetadrho_func(beta_max,0,c1,c2)*(1 + (delta - 1)*eta)
    zeta_flat = -ffmB.dbetadrho_func(beta_max,0,c1,c2)*d_tilde*(1 + (delta - 1)*eta) + ffmB.beta_func(beta_max,0,c1,c2)*d_tilde*eta*(delta - 1)

    zeta_root = np.sqrt(zeta_flat**2 - 4*d_tilde*ffmB.dbetadrho_func(beta_max,0,c1,c2)*(1 + (delta - 1)*eta))
    
    
    zeta_1 = (-eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * delta + ffmB.beta_func(beta_max,0,c1,c2) * eta * delta + eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) - ffmB.beta_func(beta_max,0,c1,c2) * eta - 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) + np.sqrt(eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 * delta ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta ** 2 + eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 * delta ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 * delta + 4 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta - 2 * eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 * delta + eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) + eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 - 4 * eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta + 4 * eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2))) / ffmB.dbetadrho_func(beta_max,0,c1,c2) / (eta * delta - eta + 1) * d_tilde / 2

    zeta_2 = -(eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * delta - ffmB.beta_func(beta_max,0,c1,c2) * eta * delta - eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) + ffmB.beta_func(beta_max,0,c1,c2) * eta + np.sqrt(eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 * delta ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta ** 2 + eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 * delta ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 * delta + 4 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta - 2 * eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 * delta + eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) ** 2 - 2 * eta ** 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) + eta ** 2 * ffmB.beta_func(beta_max,0,c1,c2) ** 2 - 4 * eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2) * delta + 4 * eta * ffmB.dbetadrho_func(beta_max,0,c1,c2) * ffmB.beta_func(beta_max,0,c1,c2)) + 2 * ffmB.dbetadrho_func(beta_max,0,c1,c2)) / ffmB.dbetadrho_func(beta_max,0,c1,c2) / (eta * delta - eta + 1) * d_tilde / 2

    if np.isnan(zeta_1) or np.isnan(zeta_2):
        zeta_1 = 0
        zeta_2 = 10

    zeta_lower = min(zeta_1, zeta_2)
    zeta_upper = max(zeta_1, zeta_2)
    temp_extinction_zetas = []
    for i in range(resolution + 1):
        rho = rho_max*i/resolution
        if rho != 1:
            extinction_condition = ((-b + d) * ffmB.beta_func(beta_max,rho,c1,c2) + b * q * (d + alpha + gamma)) * (d + alpha + gamma) / (-(b - d) * (1 + (delta - 1) * eta) * ffmB.beta_func(beta_max,rho,c1,c2) + b * q * (d + alpha + gamma)) / (-1 + rho)
        else:
            extinction_condition = (((-b + d) * ffmB.beta_func(beta_max,rho,c1,c2) + b * q * (d + alpha + gamma)) * (d + alpha + gamma))
        temp_extinction_zetas.append(extinction_condition)
    
    zeta_limit = max(temp_extinction_zetas)
#     print(zeta_limit)
    if zeta_limit < zeta_upper:
        zetas = [1*zeta_limit*i/resolution for i in range(resolution + 1)]
    else:
        zetas = [1.1*zeta_upper*i/resolution for i in range(resolution + 1)]
#     

#     print(zeta_lower, zeta_upper)
evolutionary_stabilities = []
convergence_stabilities = []
evolved_rhos = []
evolved_zetas = []
evolved_rhos_rep = []
evolved_zetas_rep = []
for zeta in zetas:
    attractors, repellers = ffmB.find_sing_strats(solB, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, zeta, seed, scenario)
    for rho in attractors:
        evolved_rhos.append(rho)
        evolved_zetas.append(zeta)
    for rho in repellers:
        evolved_rhos_rep.append(rho)
        evolved_zetas_rep.append(zeta)

prevs, prevs_no_evos = ffmB.calculate_disease_prevalance(solB, evolved_rhos, evolved_zetas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)
R0s, R0s_no_evos = ffmB.R0_model_B(solB, evolved_rhos, evolved_zetas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)

ax[0][1].plot(evolved_zetas, evolved_rhos)
# plt.scatter(evolved_zetas, simulated_rhos)
if delta == 0 and eta == 1.0:
    ax[0][1].plot([zeta_calc, zeta_calc], [0, 1.01], "k--")
    ax[0][1].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
else:
    ax[0][1].plot([zeta_lower, zeta_lower], [0, 1.01], "k--")
    ax[0][1].plot([zeta_upper, zeta_upper], [0, 1.01], "k--")
    ax[0][1].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
ax[0][1].set_ylim([0,1.01])
ax[0][1].set_xlim([0,1.01*zetas[-1]])
# ax[0][1].set_ylabel(r"Evolved false negativity, $\rho^*$")
# ax[1][1].set_xlabel(r"Testing rate, $\zeta$")
ax[0][1].set_title("Model B")
ax[0][1].text(0.03, 0.90,labels[0][1], transform=ax[0][1].transAxes)

ax[1][1].plot(evolved_zetas, prevs)
ax[1][1].plot(evolved_zetas, prevs_no_evos)
ax[1][1].set_ylim([0,0.2])
ax[1][1].set_xlim([0,1.01*zetas[-1]])
if delta == 0 and eta == 1.0:
    ax[1][1].plot([zeta_calc, zeta_calc], [0, 1.01], "k--")
    ax[1][1].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
else:
    ax[1][1].plot([zeta_lower, zeta_lower], [0, 1.01], "k--")
    ax[1][1].plot([zeta_upper, zeta_upper], [0, 1.01], "k--")
    ax[1][1].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
# ax[1][1].set_ylabel("Pathogen Prevalance")
ax[1][1].text(0.03, 0.90,labels[1][1], transform=ax[1][1].transAxes)

ax[2][1].plot(evolved_zetas, R0s)
ax[2][1].plot(evolved_zetas, R0s_no_evos)
ax[2][1].set_xlim([0,1.01*zetas[-1]])
ax[2][1].set_ylim([0,R0_upper_limit])
ax[2][1].plot([0,zetas[-1]], [1, 1], "k--")
if delta == 0 and eta == 1.0:
    ax[2][1].plot([zeta_calc, zeta_calc], [0, R0_upper_limit], "k--")
    ax[2][1].plot([zeta_limit, zeta_limit], [0, R0_upper_limit], "k-.")
else:
    ax[2][1].plot([zeta_lower, zeta_lower], [0, R0_upper_limit], "k--")
    ax[2][1].plot([zeta_upper, zeta_upper], [0, R0_upper_limit], "k--")
    ax[2][1].plot([zeta_limit, zeta_limit], [0, R0_upper_limit], "k-.")
# ax[2][1].set_ylabel(r"Basic Reproduction Number, $R_0$")
ax[2][1].set_xlabel(r"Testing rate, $\zeta$")
ax[2][1].text(0.03, 0.90,labels[2][1], transform=ax[2][1].transAxes)
# plt.savefig("../outputs/evolution_Model_A_and_B_imperfect.png", bbox_inches = "tight")
# plt.savefig("../outputs/evolution_Model_A_and_B_imperfect.pdf", bbox_inches = "tight")
# plt.close()

evolutionary_stabilities = []
convergence_stabilities = []
evolved_rhos = []
evolved_zetas = []
evolved_rhos_rep = []
evolved_zetas_rep = []
zeta_threshold_model_C = (ffmC.dbetadrho_func(beta_max,0,c1,c2)*(d + alpha + gamma))/(ffmC.dbetadrho_func(beta_max,0,c1,c2)*(-1 -delta*eta + eta) + eta*ffmC.beta_func(beta_max,0,c1,c2)*(delta - 1))
# print(zetas)
for zeta in zetas:
    attractors, repellers = ffmC.find_sing_strats(solC, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, zeta, seed, scenario)
#     print(attractors)
    for rho in attractors:
        evolved_rhos.append(rho)
        evolved_zetas.append(zeta)
    for rho in repellers:
        evolved_rhos_rep.append(rho)
        evolved_zetas_rep.append(zeta)

# flop
prevs, prevs_no_evos = ffmC.calculate_disease_prevalance(solC, evolved_rhos, evolved_zetas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)
# print(prevs)
R0s, R0s_no_evos = ffmC.R0_model_C(solC, evolved_rhos, evolved_zetas, b, d, q, alpha, gamma, c1, c2, delta, eta, seed, beta_max)

ax[0][2].plot(evolved_zetas, evolved_rhos)
ax[0][2].plot([zeta_threshold_model_C, zeta_threshold_model_C], [0, 1.01], "k--")
# plt.scatter(evolved_zetas, simulated_rhos)
# if delta == 0 and eta == 1.0:
#     ax[0][2].plot([zeta_calc, zeta_calc], [0, 1.01], "k--")
#     ax[0][2].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
# else:
#     ax[0][2].plot([zeta_lower, zeta_lower], [0, 1.01], "k--")
#     ax[0][2].plot([zeta_upper, zeta_upper], [0, 1.01], "k--")
#     ax[0][2].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
ax[0][2].set_ylim([0,1.01])
ax[0][2].set_xlim([0,1.01*zetas[-1]])
# ax[0][1].set_ylabel(r"Evolved false negativity, $\rho^*$")
# ax[1][1].set_xlabel(r"Testing rate, $\zeta$")
ax[0][2].set_title("Model C")
ax[0][2].text(0.03, 0.90,labels[0][2], transform=ax[0][2].transAxes)

ax[1][2].plot(evolved_zetas, prevs)
ax[1][2].plot(evolved_zetas, prevs_no_evos)
ax[1][2].set_ylim([0,0.2])
ax[1][2].set_xlim([0,1.01*zetas[-1]])
ax[1][2].plot([zeta_threshold_model_C, zeta_threshold_model_C], [0, 1.01], "k--")
# if delta == 0 and eta == 1.0:
#     ax[1][2].plot([zeta_calc, zeta_calc], [0, 1.01], "k--")
#     ax[1][2].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
# else:
#     ax[1][2].plot([zeta_lower, zeta_lower], [0, 1.01], "k--")
#     ax[1][2].plot([zeta_upper, zeta_upper], [0, 1.01], "k--")
#     ax[1][2].plot([zeta_limit, zeta_limit], [0, 1.01], "k-.")
# ax[1][1].set_ylabel("Pathogen Prevalance")
ax[1][2].text(0.03, 0.90,labels[1][2], transform=ax[1][2].transAxes)

ax[2][2].plot(evolved_zetas, R0s)
ax[2][2].plot(evolved_zetas, R0s_no_evos)
ax[2][2].set_xlim([0,1.01*zetas[-1]])
ax[2][2].set_ylim([0,R0_upper_limit])
ax[2][2].plot([0,zetas[-1]], [1, 1], "k--")
ax[2][2].plot([zeta_threshold_model_C, zeta_threshold_model_C], [0, R0_upper_limit], "k--")
# if delta == 0 and eta == 1.0:
#     ax[2][2].plot([zeta_calc, zeta_calc], [0, R0_upper_limit], "k--")
#     ax[2][2].plot([zeta_limit, zeta_limit], [0, R0_upper_limit], "k-.")
# else:
#     ax[2][2].plot([zeta_lower, zeta_lower], [0, R0_upper_limit], "k--")
#     ax[2][2].plot([zeta_upper, zeta_upper], [0, R0_upper_limit], "k--")
#     ax[2][2].plot([zeta_limit, zeta_limit], [0, R0_upper_limit], "k-.")
ax[2][1].set_ylabel(r"Basic Reproduction Number, $R_0$")
ax[2][2].set_xlabel(r"Testing rate, $\zeta$")
ax[2][2].text(0.03, 0.90,labels[2][2], transform=ax[2][2].transAxes)
plt.savefig("../outputs/evolution_Models_imperfect.png", bbox_inches = "tight")
plt.savefig("../outputs/evolution_Models_imperfect.pdf", bbox_inches = "tight")
plt.close()