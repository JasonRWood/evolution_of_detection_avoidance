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
    zeta,
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
            zeta,
            seed,
            S_density=4.0,
            Q_density=0.0,
            A_density=0.0,
            U_density=1.0,
        )

        fit_grad = fitness_gradient_rho(
            rho_value,
            beta_max,
            c1,
            c2,
            y,
            eta,
            zeta,
            delta,
            d,
            alpha,
            gamma,
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
        #         print(zeta, fit_grad_vec)
        if fit_grad_vec[0] == -1:
            attractors.append(rhos[0])
        if fit_grad_vec[-1] == 1:
            attractors.append(rhos[-1])

    #         print(fit_grad_vec)
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
                zeta,
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
                zeta,
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


def fitness_gradient_rho(
    rho, beta_max, c1, c2, y, eta, zeta, delta, d, alpha, gamma, scenario=1
):

    S = y[0]
    Q = y[1] * (y[1] > 1e-3)
    A = y[2] * (y[2] > 1e-3)
    U = y[3] * (y[3] > 1e-3)

    beta = beta_func(beta_max, rho, c1, c2)
    dbetadrho = dbetadrho_func(beta_max, rho, c1, c2)

    #     if delta == 0 and eta == 1:
    #         dudt = beta*S*U - zeta*(1 - rho) - (d + alpha + gamma)
    #         if dudt < 0:
    #             Q = 0
    #             A = 0
    #             U = 0
    #     print(zeta, Q, A, U)
    #     fit_grad_rho = dbetadrho*((1 - rho)*(1 - eta) + rho) - (1 - eta)*beta + beta
    #     fit_grad_rho = dbetadrho*(rho*eta - eta + 1) + beta*eta

    #     if scenario == 1:
    #         fit_grad_rho = dbetadrho*S + zeta
    #     elif scenario == 2:
    #         fit_grad_rho = ((zeta * (1 - rho) + d + alpha + gamma) * (zeta * (-1 + rho) * (-1 + eta) + d + alpha + gamma) * dbetadrho + eta * zeta * beta * (d + alpha + gamma)) * S / (d + alpha + gamma) / (zeta * (1 - rho) + d + alpha + gamma) ** 2

    #     elif scenario == 3:
    #         fit_grad_rho = ((zeta * (1 - rho) + d + alpha + gamma) * (-(1 + (delta - 1) * eta) * (-1 + rho) * zeta + d + alpha + gamma) * dbetadrho - eta * zeta * beta * (delta - 1) * (d + alpha + gamma)) / (d + alpha + gamma) * S / (zeta * (1 - rho) + d + alpha + gamma) ** 2

    #     if Q + A + U != 0:
    #         d_tilde = d + alpha + gamma
    #         denominator = (d_tilde)*((zeta*(1 - rho) + d_tilde)**2)
    fit_grad_rho = (
        (
            (-(-1 + rho) * (1 + (delta - 1) * eta) * zeta + d + alpha + gamma)
            * (zeta * (1 - rho) + d + alpha + gamma)
            * dbetadrho
            - eta * zeta * beta * (delta - 1) * (d + alpha + gamma)
        )
        * S
        / (d + alpha + gamma)
        / (zeta * (1 - rho) + d + alpha + gamma) ** 2
    )

    #         fit_grad_rho = dbetadrho*S + zeta
    #     else:
    #         print(zeta,rho)
    #         fit_grad_rho = 0
    fit_grad_rho = (fit_grad_rho > 0) * 1 + (fit_grad_rho < 0) * -1

    return fit_grad_rho


def calculate_stabilities(
    rho_t, beta_max, c1, c2, eta, delta, d, alpha, gamma, zeta, scenario=1
):
    derivative_distance = 1e-1

    rho = rho_t - derivative_distance / 2
    rho_2 = rho_t + derivative_distance / 2

    beta = beta_func(beta_max, rho, c1, c2)
    beta_2 = beta_func(beta_max, rho_2, c1, c2)

    dbetadrho = dbetadrho_func(beta_max, rho, c1, c2)
    dbetadrho_2 = dbetadrho_func(beta_max, rho_2, c1, c2)

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
        zeta,
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
        zeta,
        seed,
        S_density=4.0,
        Q_density=0.0,
        A_density=0.0,
        U_density=1.0,
    )

    S = y[0]
    S_2 = y_2[0]

    #     if scenario == 1:
    #         fit_grad_rho = dbetadrho*S + zeta
    #     elif scenario == 2:
    #         fit_grad_rho = ((zeta * (1 - rho) + d + alpha + gamma) * (zeta * (-1 + rho) * (-1 + eta) + d + alpha + gamma) * dbetadrho + eta * zeta * beta * (d + alpha + gamma)) * S / (d + alpha + gamma) / (zeta * (1 - rho) + d + alpha + gamma) ** 2

    #     elif scenario == 3:
    #         fit_grad_rho = ((zeta * (1 - rho) + d + alpha + gamma) * (-(1 + (delta - 1) * eta) * (-1 + rho) * zeta + d + alpha + gamma) * dbetadrho - eta * zeta * beta * (delta - 1) * (d + alpha + gamma)) / (d + alpha + gamma) * S / (zeta * (1 - rho) + d + alpha + gamma) ** 2
    denominator = (d + alpha + gamma) * ((zeta * (1 - rho) + d + alpha + gamma) ** 2)
    fit_grad_rho = (
        (
            (zeta * (1 - rho) + d + alpha + gamma)
            * (-(1 + (delta - 1) * eta) * (rho - 1) * zeta + d + alpha + gamma)
            * dbetadrho
            - eta * zeta * beta * (zeta - 1) * (d + alpha + gamma)
        )
        * S
    ) / denominator
    #     fit_grad_rho = dbetadrho*S + zeta
    #     if scenario == 1:
    #         fit_grad_rho_res = dbetadrho*S_2 + zeta
    #     elif scenario == 2:
    #         fit_grad_rho_res = ((zeta * (1 - rho) + d + alpha + gamma) * (zeta * (-1 + rho) * (-1 + eta) + d + alpha + gamma) * dbetadrho + eta * zeta * beta * (d + alpha + gamma)) * S_2 / (d + alpha + gamma) / (zeta * (1 - rho) + d + alpha + gamma) ** 2

    #     elif scenario == 3:
    #         fit_grad_rho_res = ((zeta * (1 - rho) + d + alpha + gamma) * (-(1 + (delta - 1) * eta) * (-1 + rho) * zeta + d + alpha + gamma) * dbetadrho - eta * zeta * beta * (delta - 1) * (d + alpha + gamma)) / (d + alpha + gamma) * S_2 / (zeta * (1 - rho) + d + alpha + gamma) ** 2

    denominator = (d + alpha + gamma) * ((zeta * (1 - rho) + d + alpha + gamma) ** 2)
    fit_grad_rho_res = (
        (
            (zeta * (1 - rho) + d + alpha + gamma)
            * (-(1 + (delta - 1) * eta) * (rho - 1) * zeta + d + alpha + gamma)
            * dbetadrho
            - eta * zeta * beta * (zeta - 1) * (d + alpha + gamma)
        )
        * S_2
    ) / denominator
    #     fit_grad_rho_res = dbetadrho*S_2 + zeta
    #     if scenario == 1:
    #         fit_grad_rho_mut = dbetadrho_2*S + zeta
    #     elif scenario == 2:
    #         fit_grad_rho_mut = ((zeta * (1 - rho_2) + d + alpha + gamma) * (zeta * (-1 + rho_2) * (-1 + eta) + d + alpha + gamma) * dbetadrho_2 + eta * zeta * beta_2 * (d + alpha + gamma)) * S / (d + alpha + gamma) / (zeta * (1 - rho_2) + d + alpha + gamma) ** 2

    #     elif scenario == 3:
    #         fit_grad_rho_mut = ((zeta * (1 - rho_2) + d + alpha + gamma) * (-(1 + (delta - 1) * eta) * (-1 + rho_2) * zeta + d + alpha + gamma) * dbetadrho_2 - eta * zeta * beta_2 * (delta - 1) * (d + alpha + gamma)) / (d + alpha + gamma) * S / (zeta * (1 - rho_2) + d + alpha + gamma) ** 2

    denominator = (d + alpha + gamma) * ((zeta * (1 - rho_2) + d + alpha + gamma) ** 2)
    fit_grad_rho_mut = (
        (
            (zeta * (1 - rho_2) + d + alpha + gamma)
            * (-(1 + (delta - 1) * eta) * (rho_2 - 1) * zeta + d + alpha + gamma)
            * dbetadrho_2
            - eta * zeta * beta_2 * (zeta - 1) * (d + alpha + gamma)
        )
        * S
    ) / denominator
    #     fit_grad_rho_mut = dbetadrho_2*S + zeta
    E = (fit_grad_rho_mut - fit_grad_rho) / derivative_distance
    M = (fit_grad_rho_res - fit_grad_rho) / derivative_distance

    return (E, E + M)


def calculate_disease_prevalance(
    sol,
    evolved_rhos,
    evolved_zetas,
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
        zeta = evolved_zetas[i]
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
            zeta,
            seed,
            S_density=1.0,
            Q_density=1.0,
            A_density=1.0,
            U_density=1.0,
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
            zeta,
            seed,
            S_density=1.0,
            Q_density=1.0,
            A_density=1.0,
            U_density=1.0,
        )
        total_pop = sum(y)
        disease_pop = total_pop - y[0] - y[-1]
        #         print(y)
        temp_prevalance = disease_pop / total_pop
        prevalances.append(temp_prevalance)
        prevalances_no_evo.append((sum(y2) - y2[0] - y2[-1]) / sum(y2))
    return prevalances, prevalances_no_evo


def R0_model_B(
    sol,
    evolved_rhos,
    evolved_zetas,
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
    R0s_no_evos = []
    for i, rho in enumerate(evolved_rhos):
        beta = beta_func(beta_max, rho, c1, c2)
        zeta = evolved_zetas[i]

        R0 = beta_func(beta_max, rho, c1, c2) * (b - d) / b / q * (-(-1 + rho) * (1 + (delta - 1) * eta) * zeta + d + gamma + alpha) / (d + alpha + gamma) / (zeta * (1 - rho) + d + alpha + gamma)
        
        R0_no_evo = beta_func(beta_max, 0, c1, c2) * (b - d) / b / q * (-(-1 + 0) * (1 + (delta - 1) * eta) * zeta + d + gamma + alpha) / (d + alpha + gamma) / (zeta * (1 - 0) + d + alpha + gamma)
        
#         R0 = (
#             (
#                 math.sqrt(
#                     (b - d) ** 2 * beta_func(beta_max, rho, c1, c2) ** 2
#                     - 4
#                     * (b - d)
#                     * q
#                     * (delta * eta - eta + 0.1e1 / 0.2e1)
#                     * b
#                     * (-1 + rho)
#                     * zeta
#                     * beta_func(beta_max, rho, c1, c2)
#                     + zeta ** 2 * b ** 2 * q ** 2 * (-1 + rho) ** 2
#                 )
#                 + (b - d) * beta_func(beta_max, rho, c1, c2)
#                 + ((-1 + rho) * zeta - 2 * d - 2 * gamma - 2 * alpha) * q * b
#             )
#             / b
#             / q
#             / 2
#         ) + 1

#         R0_no_evo = (
#             (
#                 math.sqrt(
#                     (b - d) ** 2 * beta_func(beta_max, 0, c1, c2) ** 2
#                     - 4
#                     * (b - d)
#                     * q
#                     * (delta * eta - eta + 0.1e1 / 0.2e1)
#                     * b
#                     * (-1 + 0)
#                     * zeta
#                     * beta_func(beta_max, 0, c1, c2)
#                     + zeta ** 2 * b ** 2 * q ** 2 * (-1 + 0) ** 2
#                 )
#                 + (b - d) * beta_func(beta_max, 0, c1, c2)
#                 + ((-1 + 0) * zeta - 2 * d - 2 * gamma - 2 * alpha) * q * b
#             )
#             / b
#             / q
#             / 2
#         ) + 1
        R0s.append(R0)
        R0s_no_evos.append(R0_no_evo)
    return R0s, R0s_no_evos
