def calculate_E(rho, beta, dbetadrho_2, delta, eta, sigma):

    E_model_A = (
        dbetadrho_2 / beta
        - 2
        * (delta - 1) ** 2
        * eta ** 2
        * sigma ** 2
        / (-1 + eta * (delta - 1) * (-1 + rho) * sigma) ** 2
    )

    return E_model_A
