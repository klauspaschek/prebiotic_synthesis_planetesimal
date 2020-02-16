import numpy as np
import scipy.optimize



####
# Calculate approximate boiling point temperature T of water as function of
# pressure P by using an approximation for the pressure temperature dependence
# P(T) = f(T). By rearranging the formula as P(T) - f(T) = g(T) = 0 and finding
# the root using the Newton-Raphson method this formula is solved for T.
####
# Parameters: pressure    : float
#                 Pressure in unit [bar].
####
# Returns:    temperature : float
#                 Temperature of boiling point of water in Kelvin.
def boiling_point(pressure):
    # Parameters of approximating function for temperature pressure dependence
    a = -6094.4642
    b = 21.1249952
    c = -2.7245552e-2
    d = 1.6853396e-5
    e = 2.4575506

    # Rearranged function P(T) - f(T)
    def f(T):
        return pressure - np.exp(a / T + b + c * T + d * T**2 + \
                                 e * np.log(T)) * 1e-5

    # Derivative of f(T)
    def fprime(T):
        return -np.exp(a / T + b + c * T + d * T**2 + e * np.log(T)) * \
               (-a / T**2 + c + 2 * d * T + e / T) * 1e-5

    temperature = scipy.optimize.newton(f, 1000, fprime = fprime)

    return temperature
