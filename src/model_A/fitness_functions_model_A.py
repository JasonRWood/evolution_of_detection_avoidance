import numpy as np
import math


def find_sing_strats(
    sol,
    depth,
    depth_max,
    resolution,
    rho_lower,
    rho_upper,
    beta_max,
    c1,
    c2,
    b,
    q,
    d,
    alpha,
    gamma,
    delta,
    eta,
    sigma,
    seed,
    scenario,
):

    attractors = []
    repellers = []
    poss_attractors = []
    poss_repellers = []
    fit_grad_vec = []
    rhos = []
    resolution_inner = resolution * (1 + 9 * (depth == 0))
    for i in range(resolution_inner + 1):
        rho_value = rho_lower + (rho_upper - rho_lower) * i / (resolution_inner)
        #         print(rho_value)
        rhos.append(rho_value)
        beta_value = beta_func(beta_max, rho_value, c1, c2)

        y = sol.eco_dynamics(
            beta_value,
            rho_value,
            b,
            q,
            d,
            alpha,
            gamma,
            c1,
            c2,
            delta,
            eta,
            sigma,
            seed,
            S_density=4.0,
            Q_density=1.0,
            A_density=1.0,
            U_density=1.0,
        )

        fit_grad = fitness_gradient_rho(
            rho_value,
            beta_max,
            c1,
            c2,
            y[0],
            eta,
            delta,
            d,
            alpha,
            gamma,
            sigma,
            scenario=scenario,
        )
        fit_grad_vec.append(fit_grad)

        if i != 0:
            if fit_grad_vec[-2] == 1 and (
                fit_grad_vec[-1] == -1 or fit_grad_vec[-1] == 0
            ):
                poss_attractors.append([rhos[-2], rhos[-1]])
            if fit_grad_vec[-2] == -1 and (
                fit_grad_vec[-1] == 1 or fit_grad_vec[-1] == 0
            ):
                poss_repellers.append([rhos[-2], rhos[-1]])

    #     print(depth, poss_attractors, fit_grad_vec)
    if depth == 0:
        #         if c1 == 0.45 and c2 == 0.0:
        #             print(fit_grad_vec)
        if fit_grad_vec[0] == -1:
            attractors.append(rhos[0])
        if fit_grad_vec[-1] == 1:
            attractors.append(rhos[-1])

    if depth != depth_max:
        for pair in poss_attractors:
            temp_attractors, temp_repellers = find_sing_strats(
                sol,
                depth + 1,
                depth_max,
                resolution,
                pair[0],
                pair[1],
                beta_max,
                c1,
                c2,
                b,
                q,
                d,
                alpha,
                gamma,
                delta,
                eta,
                sigma,
                seed,
                scenario,
            )
            for val in temp_attractors:
                attractors.append(val)
            for val in temp_repellers:
                repellers.append(val)

        for pair in poss_repellers:
            temp_attractors, temp_repellers = find_sing_strats(
                sol,
                depth + 1,
                depth_max,
                resolution,
                pair[0],
                pair[1],
                beta_max,
                c1,
                c2,
                b,
                q,
                d,
                alpha,
                gamma,
                delta,
                eta,
                sigma,
                seed,
                scenario,
            )
            for val in temp_attractors:
                attractors.append(val)
            for val in temp_repellers:
                repellers.append(val)

    else:
        for pair in poss_attractors:
            attractors.append((pair[0] + pair[1]) / 2)
        for pair in poss_repellers:
            repellers.append((pair[0] + pair[1]) / 2)

    return attractors, repellers


def beta_func(beta_max, rho, c1, c2):

    if c2 == 0.0:
        beta = beta_max * ((1 - c1) + c1 * (1 - rho))

    else:
        beta = beta_max * (
            (1 - c1) + c1 * ((1 - np.exp((1 - rho) * c2)) / (1 - np.exp(c2)))
        )

    return beta


def dbetadrho_func(beta_max, rho, c1, c2):

    if c2 == 0.0:
        dbetadrho = -beta_max * c1
    else:
        dbetadrho = beta_max * (c1 * c2 * np.exp((1 - rho) * c2)) / (1 - np.exp(c2))
    return dbetadrho


def d2betadrho2_func(beta_max, rho, c1, c2):

    if c2 == 0.0:

        d2betadrho2 = 0

    else:

        d2betadrho2 = (
            -beta_max * (c1 * c2 ** 2 * np.exp((1 - rho) * c2)) / (1 - np.exp(c2))
        )

    return d2betadrho2


def fitness_gradient_rho(
    rho, beta_max, c1, c2, S, eta, delta, d, alpha, gamma, sigma, scenario=1
):

    beta = beta_func(beta_max, rho, c1, c2)
    dbetadrho = dbetadrho_func(beta_max, rho, c1, c2)

    #     fit_grad_rho = dbetadrho*((1 - rho)*(1 - eta) + rho) - (1 - eta)*beta + beta
    #     fit_grad_rho = dbetadrho*(rho*eta - eta + 1) + beta*eta

    if scenario == 1:
        fit_grad_rho = sigma * beta * S + (sigma * rho + (1 - sigma)) * dbetadrho * S
    elif scenario == 2:
        fit_grad_rho = (
            S
            * ((1 + eta * (-1 + rho) * sigma) * dbetadrho + beta * eta * sigma)
            / (d + alpha + gamma)
        )
    elif scenario == 3:
        fit_grad_rho = (
            -S
            * (
                (-1 + eta * (delta - 1) * (-1 + rho) * sigma) * dbetadrho
                + beta * eta * sigma * (delta - 1)
            )
            / (d + alpha + gamma)
        )

    fit_grad_rho = (fit_grad_rho > 0) * 1 + (fit_grad_rho < 0) * -1

    return fit_grad_rho


def calculate_stabilities(
    rho_t, beta_max, c1, c2, eta, delta, d, alpha, gamma, sigma, scenario=1
):
    derivative_distance = 1e-1

    rho = rho_t - derivative_distance / 2
    rho_2 = rho_t + derivative_distance / 2

    beta_value = beta_func(beta_max, rho, c1, c2)
    beta_2 = beta_func(beta_max, rho_2, c1, c2)

    dbetadrho = dbetadrho_func(beta_max, rho, c1, c2)
    dbetadrho_2 = dbetadrho_func(beta_max, rho_2, c1, c2)

    y = sol.eco_dynamics(
        beta_value,
        rho,
        b,
        q,
        d,
        alpha,
        gamma,
        c1,
        c2,
        delta,
        eta,
        sigma,
        seed,
        S_density=4.0,
        Q_density=1.0,
        A_density=1.0,
        U_density=1.0,
    )

    y_2 = sol.eco_dynamics(
        beta_2,
        rho_2,
        b,
        q,
        d,
        alpha,
        gamma,
        c1,
        c2,
        delta,
        eta,
        sigma,
        seed,
        S_density=4.0,
        Q_density=1.0,
        A_density=1.0,
        U_density=1.0,
    )

    S = y[0]
    S_2 = y_2[0]

    if scenario == 1:
        fit_grad_rho = (
            sigma * beta_value * S + (sigma * rho + (1 - sigma)) * dbetadrho * S
        )
    elif scenario == 2:
        fit_grad_rho = (
            S
            * ((1 + eta * (-1 + rho) * sigma) * dbetadrho + beta_value * eta * sigma)
            / (d + alpha + gamma)
        )
    elif scenario == 3:
        fit_grad_rho = (
            -S
            * (
                (-1 + eta * (delta - 1) * (-1 + rho) * sigma) * dbetadrho
                + beta_value * eta * sigma * (delta - 1)
            )
            / (d + alpha + gamma)
        )

    if scenario == 1:
        fit_grad_rho_res = (
            sigma * beta_value * S_2 + (sigma * rho + (1 - sigma)) * dbetadrho * S_2
        )
    elif scenario == 2:
        fit_grad_rho_res = (
            S_2
            * ((1 + eta * (-1 + rho) * sigma) * dbetadrho + beta_value * eta * sigma)
            / (d + alpha + gamma)
        )
    elif scenario == 3:
        fit_grad_rho_res = (
            -S_2
            * (
                (-1 + eta * (delta - 1) * (-1 + rho) * sigma) * dbetadrho
                + beta_value * eta * sigma * (delta - 1)
            )
            / (d + alpha + gamma)
        )

    if scenario == 1:
        fit_grad_rho_mut = (
            sigma * beta_2 * S + (sigma * rho_2 + (1 - sigma)) * dbetadrho_2 * S
        )
    elif scenario == 2:
        fit_grad_rho_mut = (
            S
            * ((1 + eta * (-1 + rho_2) * sigma) * dbetadrho_2 + beta_2 * eta * sigma)
            / (d + alpha + gamma)
        )
    elif scenario == 3:
        fit_grad_rho_mut = (
            -S
            * (
                (-1 + eta * (delta - 1) * (-1 + rho_2) * sigma) * dbetadrho_2
                + beta_2 * eta * sigma * (delta - 1)
            )
            / (d + alpha + gamma)
        )

    E = (fit_grad_rho_mut - fit_grad_rho) / derivative_distance
    M = (fit_grad_rho_res - fit_grad_rho) / derivative_distance

    return (E, E + M)


def calculate_disease_prevalance(
    sol,
    evolved_rhos,
    evolved_sigmas,
    b,
    d,
    q,
    alpha,
    gamma,
    c1,
    c2,
    delta,
    eta,
    seed,
    beta_max,
):

    prevalances = []
    prevalances_no_evo = []

    for i, rho in enumerate(evolved_rhos):
        beta = beta_func(beta_max, rho, c1, c2)
        sigma = evolved_sigmas[i]
        y = sol.eco_dynamics(
            beta,
            rho,
            b,
            q,
            d,
            alpha,
            gamma,
            c1,
            c2,
            delta,
            eta,
            sigma,
            seed,
            S_density=1.0,
            Q_density= 0.10,
            A_density= 0.10,
            U_density= 0.10,
        )

        y2 = sol.eco_dynamics(
            beta_func(beta_max, 0, c1, c2),
            0,
            b,
            q,
            d,
            alpha,
            gamma,
            c1,
            c2,
            delta,
            eta,
            sigma,
            seed,
            S_density=1.0,
            Q_density= 0.10,
            A_density= 0.10,
            U_density= 0.10,
        )
        total_pop = sum(y)
        disease_pop = total_pop - y[0] - y[-1]
        #         print(y)
        temp_prevalance = disease_pop / total_pop
        prevalances.append(temp_prevalance)
        
        total_pop = sum(y2)
        disease_pop = total_pop - y2[0] - y2[-1]
        #         print(y)
        temp_prevalance = disease_pop / total_pop
        prevalances_no_evo.append(temp_prevalance)
    return prevalances, prevalances_no_evo


def R0s_model_A(
    sol,
    evolved_rhos,
    evolved_sigmas,
    b,
    d,
    q,
    alpha,
    gamma,
    c1,
    c2,
    delta,
    eta,
    seed,
    beta_max,
):

    R0s = []
    R0s_no_evo = []
    for i, rho in enumerate(evolved_rhos):
        beta = beta_func(beta_max, rho, c1, c2)
        sigma = evolved_sigmas[i]

        R0 = -(b - d) / b / q * (-1 + eta * (delta - 1) * (-1 + rho) * sigma) * beta_func(beta_max, rho, c1, c2) / (d + alpha + gamma)
        
        R0_no_evo = -(b - d) / b / q * (-1 + eta * (delta - 1) * (-1 + 0) * sigma) * beta_func(beta_max, 0, c1, c2) / (d + alpha + gamma)
#         R0 = (
#             -(
#                 beta_func(beta_max, rho, c1, c2) * b * delta * eta * rho * sigma
#                 - beta_func(beta_max, rho, c1, c2) * d * delta * eta * rho * sigma
#                 - beta_func(beta_max, rho, c1, c2) * b * delta * eta * sigma
#                 - beta_func(beta_max, rho, c1, c2) * b * eta * rho * sigma
#                 + beta_func(beta_max, rho, c1, c2) * d * delta * eta * sigma
#                 + beta_func(beta_max, rho, c1, c2) * d * eta * rho * sigma
#                 + beta_func(beta_max, rho, c1, c2) * b * eta * sigma
#                 - beta_func(beta_max, rho, c1, c2) * d * eta * sigma
#                 + alpha * b * q
#                 + d * b * q
#                 + gamma * b * q
#                 - beta_func(beta_max, rho, c1, c2) * b
#                 + beta_func(beta_max, rho, c1, c2) * d
#             )
#             / b
#             / q
#         ) + 1

#         R0_no_evo = (
#             -(
#                 beta_func(beta_max, 0, c1, c2) * b * delta * eta * 0 * sigma
#                 - beta_func(beta_max, 0, c1, c2) * d * delta * eta * 0 * sigma
#                 - beta_func(beta_max, 0, c1, c2) * b * delta * eta * sigma
#                 - beta_func(beta_max, 0, c1, c2) * b * eta * 0 * sigma
#                 + beta_func(beta_max, 0, c1, c2) * d * delta * eta * sigma
#                 + beta_func(beta_max, 0, c1, c2) * d * eta * 0 * sigma
#                 + beta_func(beta_max, 0, c1, c2) * b * eta * sigma
#                 - beta_func(beta_max, 0, c1, c2) * d * eta * sigma
#                 + alpha * b * q
#                 + d * b * q
#                 + gamma * b * q
#                 - beta_func(beta_max, 0, c1, c2) * b
#                 + beta_func(beta_max, 0, c1, c2) * d
#             )
#             / b
#             / q
#         ) + 1

        R0s.append(R0)
        R0s_no_evo.append(R0_no_evo)
    return R0s, R0s_no_evo
