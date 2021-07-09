"""
Define all the physical and mathematical constants that we will use
"""
global PI, HBARC, M_PROTON, M_NEUTRON, M_NUCLEON, gA, gpipn, Fpi, mpi, LambdaA, LambdaChi, c6_hat, r2A, gpinn
PI = 3.14159265359
HBARC = 197.3269718  # reduced Planck constant * speed of light in MeV * fm
M_PROTON = 938.2720813  # proton mass in MeV/c^2
M_NEUTRON = 939.5654133  # neutron mass in MeV/c^2
M_NUCLEON = (M_PROTON+M_NEUTRON) / 2.0  # avarage nucleon mass in MeV/c^2


"""
from N3LO(EM), NNLOsat.  Physics Reports503,1(2011)
"""
gA = 1.29
gpinn = 13.1066 # gpinn^2/(4*math.pi) = 13.67


"""
Values of constants from J. Men\'edez, D. Gazit and A. Schwenk, PRD86,103511(2012)
"""
gpipn = 13.05
Fpi = 92.4  # in MeV
mpi = 138.04  # in MeV
LambdaA = 1040.0  # in MeV
LambdaChi = 700  # in MeV
c6_hat = 5.83  # from V. Bernard, N. Kaiser, and U.-G. Meißner, Nucl.Phys.A615,483(1997)


"""
From P. Klos, et al., PRD88,083516(2013)
"""
# gA = 1.2670


"""
From M. Hoferichter, et al., PRD102,074018(2020)
"""
# gA = 1.27641
# Fpi = 92.28  # in MeV
# gpinn = 13.1209 # gpinn^2/(4*math.pi) = 13.7(2)

r2A = 0.47 # in fm^2
