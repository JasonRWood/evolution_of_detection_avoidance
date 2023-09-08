def calculate_E(rho, beta, dbetadrho_2, d, alpha, gamma, eta, delta, zeta):

    E_model_C = (
        (
            -2 * (delta - 1) ** 2 * eta ** 2 * zeta ** 2 * beta
            + ((1 - (delta - 1) * (-1 + rho) * eta) * zeta + d + gam + alpha) ** 2
            * dbetadrho_2
        )
        / beta
        / ((1 - (delta - 1) * (-1 + rho) * eta) * zeta + d + gam + alpha) ** 2
    )

    return E_model_B
