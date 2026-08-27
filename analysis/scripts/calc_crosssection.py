import numpy as np
from uncertainties import ufloat

c = 299792458 # [m/s]

M_LiBr = 6.941 + 79.904 # g/mol
N_A = 6.02214076 * 10**23 # mol^-1
natural_abundance_6Li = 0.075
LiBr_ratio = 31.0 # g/L
N_6Li = LiBr_ratio / M_LiBr * N_A * natural_abundance_6Li / 0.001 # [m^-3]
print(N_6Li)
tau = ufloat(143.8, 4.4) * 10**-6 # [s]

K_n = 0.025 # [eV]
M_n = 939.56542052 # MeV/c2 
beta = np.sqrt(1.0 - 1.0 / (1.0 + K_n*10**-6/M_n)**2)
v_n = beta * c # [m/s]

sigma = 1/tau / N_6Li / v_n # [m^2]
sigma = sigma * 10**28 # [barn]
print(sigma)

# ------------------------

rho_PC = 0.8761 # [g/cm^3]
M_PC = 120.19 # [g/mol]
N_PC = rho_PC / M_PC * N_A # [cm^-3]
N_proton = N_PC * 12
sigma_np = 0.33 # [barn]
sigma_np = sigma_np * 10**-24 # [cm^2]

v_n_cm = v_n * 100.0

tau_np = 1.0 / (sigma_np * N_proton * v_n_cm)
print(tau_np)
tau_corr = (1.0/tau - 1.0/tau_np)**-1
print(tau_corr)
sigma = 1/tau_corr / N_6Li / v_n
sigma = sigma * 1e28
print(sigma)

# ------------------------

sigma_n6Li = 941.0 # [barn]
sigma_n6Li = sigma_n6Li * 10**-28 # [m^2]
tau_n6Li = 1.0/(sigma_n6Li * N_6Li * v_n)
print(tau_n6Li)

tau_residual = (1.0/tau - 1.0/tau_n6Li - 1.0/tau_np)**-1
print(1.0 / tau_residual)