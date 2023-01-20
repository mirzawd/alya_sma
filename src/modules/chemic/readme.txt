- n PDE's
- m ODE's
- n+m coupled equations
- A+B ->(k+) <-(k-) C
  k+ = 4*pi*r_AB*(k_A+k_B) (k=diffusion coefficient)
  k- = rho_materiau * k+ * exp(-E/kT)

  r_AB: analytical radius = r_A + r_B = sum of geometrical radii
  r_A = [ 3/4*N_A/(pi*rho_materiau) ], N_A: number of specy atoms

Input
-----
- n diffusion coefficients k_i

IRRADIATED_MATERIAL
- (n+m) Number of atoms N_i
- (n+m) x (n+m) reaction matrix
- rho_materiau: material density (1 constant)
- Binding Energy E
- Temperature T
- Boltzmann constant = 1.380?6504 x 10-23 J/K 
                     = 8.617?343  x 10-5 eV/K 

To compute diffusion k=k0*exp(-Em/kT), for each specy:
- Em= migration energy
- k0= pre-exponential factor=1.37 x 10^-3 cm^2/s
