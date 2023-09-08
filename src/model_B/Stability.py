def calculate_E(rho, beta, dbetadrho_2, d, alpha, gamma, eta, delta, zeta):

    E_model_B = (
        1
        / ((1 - rho) * delta + d + gamma + alpha) ** 2
        * (-zeta * rho + zeta + alpha + d + gamma)
        * (
            -2
            * eta
            * delta ** 2
            * (-1 + delta)
            * (delta * eta - eta + 1)
            * (d + alpha + gamma)
            * beta
            + ((1 - rho) * delta + d + gamma + alpha)
            * dbetadrho_2
            * (
                -delta ** 2 * (-1 + rho) * eta
                + delta * (-1 + rho) * (-1 + eta)
                + d
                + gamma
                + alpha
            )
            ** 2
        )
        / beta
        / (
            -delta ** 2 * (-1 + rho) * eta
            + delta * (-1 + rho) * (-1 + eta)
            + d
            + gamma
            + alpha
        )
        / (
            -eta * zeta * (-1 + rho) * delta
            + zeta * (-1 + rho) * eta
            - zeta * rho
            + d
            + gamma
            + zeta
            + alpha
        )
    )

    return E_model_B
