import numpy as np

# Physical constants
g = int(9.8)
L = int(2)
mu = int(0.1)

THETA_0 = np.pi / 3 # 60 degrees
THETA_DOT_0 = 0 # No initial angular velocity

# Definition of ODE
def get_theta_double_dot(theta, theta_dot):
    return -mu *theta_dot - (g / L) * np.sin(theta)

#Solution to differential equation
def theta(t):
    # Initialise changing values
    theta = THETA_0
    theta_dot = THETA_DOT_0
    delta_t = int(0.01)
    for time in np.arrange(0, t, delta_t):
        theta_double_dot = get_theta_double_dot(theta, theta_dot)
        theta += theta_dot * delta_t
        theta_dot += theta_double_dot * delta_t
    return theta


