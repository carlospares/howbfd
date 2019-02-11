def wenorec(order, args):
    """ Wrapper for WENO reconstructions in module wenorec
        Input:
            order: order of reconstruction
            args: values at stencil, left to right (in array of length = order)
        Output: 
            double with numeric reconstruction
    """
    if order == 3:
        return weno3_rec(*args)
    elif order == 5:
        return weno5_rec(*args)


def weno3_rec(phim1, phi0, phip1):
    """ Performs left-biased 3rd order reconstruction
    Input:
        phim1: cell avg at I_{i-1}
        phi0: cell avg at I_i
        phip1: cell avg at I_{i+1}
    Output:
        reconstruction at x_{i+1/2}
    Remark:
        By symmetry, weno3_rec(phip1, phi0, phim1) will perform a right-biased
        reconstruction at x_{i-1/2}
    """
    eps = 1e-10
    p0 = -0.5*phim1 + 1.5*phi0
    p1 = 0.5*phi0 + 0.5*phip1
     
    beta0 = (phi0 - phim1)*(phi0 - phim1)
    beta1 = (phip1 - phi0)*(phip1 - phi0)
     
    alpha0 = (1.0/3.0) / (eps + beta0) / (eps + beta0)
    alpha1 = (2.0/3.0) / (eps + beta1) / (eps + beta1)
     
    alpha_sum = alpha0 + alpha1
     
    w0 = alpha0 / alpha_sum
    w1 = alpha1 / alpha_sum
    
    return w0 * p0 + w1 * p1


def weno5_rec(phim2, phim1, phi0, phip1, phip2):
    """ Performs left-biased 5th order reconstruction
        Cf. weno3_rec
    """
    eps = 1e-10
    p0 = (1.0/3.0) * phim2  - (7.0/6.0)*phim1 + (11.0/6.0)*phi0
    p1 = (-1.0/6.0) * phim1 + (5.0/6.0)*phi0 + (1.0/3.0)*phip1
    p2 = (1.0/3.0) * phi0 + (5.0/6.0)*phip1 - (1.0/6.0)*phip2

    beta2 = 13.0/12.0 * (phi0 - 2.0 * phip1 + phip2)*(phi0 - 2.0 * phip1 + phip2) + 0.25 * (3.0 * phi0 - 4.0 * phip1 + phip2)*(3.0 * phi0 - 4.0 * phip1 + phip2)
    beta1 = 13.0/12.0 * (phim1 - 2.0 * phi0 + phip1)*(phim1 - 2.0 * phi0 + phip1) + 0.25 * (phim1 - phip1)*(phim1 - phip1)
    beta0 = 13.0/12.0 * (phim2 - 2.0 * phim1 + phi0)*(phim2 - 2.0 * phim1 + phi0) + 0.25 * (phim2 - 4.0 * phim1 + 3.0 * phi0)*(phim2 - 4.0 * phim1 + 3.0 * phi0)

    alpha0 = 0.1 /(beta0 + eps)/(beta0 + eps)
    alpha1 = 0.6 /(beta1 + eps)/(beta1 + eps)
    alpha2 = 0.3 /(beta2 + eps)/(beta2 + eps)

    alpha_sum = alpha0 + alpha1 + alpha2

    w0 = alpha0/alpha_sum
    w1 = alpha1/alpha_sum
    w2 = alpha2/alpha_sum


    phi_half = w0 * p0 + w1 * p1 + w2 * p2

    return phi_half