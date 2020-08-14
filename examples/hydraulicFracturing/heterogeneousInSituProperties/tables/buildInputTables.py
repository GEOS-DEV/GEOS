
import numpy as np
from scipy.signal import medfilt
import matplotlib.pyplot as plt


def random_normal_fractal(z, fdim, sharpen=0):
  # Generate the values on a larger grid (N needs to be even)
  Na = len(z)
  dz = abs(z[1] - z[0])
  Nb = 4 * Na
  dim = 1.0

  # Setup wavenumber grid
  fn = 0.5 / dz
  k = np.linspace(-1.0, 1.0, Nb) * fn

  # Generate filter
  k2 = k**2
  k2 = np.maximum(k2, min(k2[k2.nonzero()]))
  F = k2**(-0.25*dim - 0.5*fdim)

  # Generate random input and apply the filter
  fV = np.fft.fft(np.random.randn(Nb)) * np.fft.fftshift(F)
  V = np.real(np.fft.ifft(fV))

  # Scale the results
  dist = (V - V.mean()) / V.std()

  # Filter values
  if (sharpen > 0):
    dist = medfilt(dist, sharpen)

  # Choose a subset of the distribution to return
  # Try centering around a local minimum
  mid = Nb//2
  tmp_dist = dist[mid-Na:mid+Na]
  offset = np.argmin(tmp_dist) - Na
  Nstart = Nb//2 - Na//2 + offset
  dist_new = dist[Nstart:Nstart+Na]

  return dist_new


def convert_to_bulk_shear_modulus(E, nu):
  K = E / (3.0 * (1 - 2.0 * nu))
  G = E / (2.0 * (1 + nu))
  return K, G


def build_geologic_model(z, z_o, density, E_mean, E_std, nu_mean, nu_std, sigma_tectonic, fdim):
  # Randomly generate in-situ E, nu
  E = E_mean + E_std * random_normal_fractal(z, fdim, sharpen=3)
  nu = nu_mean + nu_std * random_normal_fractal(z, fdim, sharpen=3)
  bulk_modulus, shear_modulus = convert_to_bulk_shear_modulus(E, nu)
  k = nu / (1.0 - nu)

  # Generate model
  model = {}
  model['x'] = [0.0]
  model['y'] = [0.0]
  model['z'] = z
  model['sigma_zz'] = density * -9.81 * (z - z_o)
  model['sigma_xx'] = model['sigma_zz'] * k
  model['sigma_yy'] = model['sigma_zz'] * k
  model['bulkModulus'] = bulk_modulus
  model['shearModulus'] = shear_modulus
  model['porePressure'] = 1000.0 * -9.81 * (z - z_o)
  # model['nu'] = nu

  # Add in tectonic stress, flip sign
  model['sigma_xx'] = -1.0 * model['sigma_xx'] - sigma_tectonic[0]
  model['sigma_yy'] = -1.0 * model['sigma_yy'] - sigma_tectonic[1]
  model['sigma_zz'] = -1.0 * model['sigma_zz'] - sigma_tectonic[2]

  return model


def write_table_files(tables, output_path, output_format='%1.5e'):
  for k in tables.keys():
    tmp = np.reshape(tables[k], (-1), order='F')
    np.savetxt('%s/%s.csv' % (output_path, k), tmp, fmt=output_format, delimiter=',')


# Config
# Output
output_path = './'

# Stress model
z = np.linspace(-150, 150, 76) + 2.0
z_o = 500.0
E_mean = 2e10
E_std = 1e9
nu_mean = 0.25
nu_std = 0.02
density = 2650.0
sigma_tectonic = [4.0e6, 3.0e6, 0.0]
fdim = 0.0
seed = 1234567

# Pumping schedule
q = 0.05
t_init = 60.0
t_ramp = 60.0
t_pump = 60 * 60.0


# Build the geologic model
np.random.seed(seed)
model = build_geologic_model(z, z_o, density, E_mean, E_std, nu_mean, nu_std, sigma_tectonic, fdim)

# Build the pumping schedule
model['flowRate'] = [0.0, 0.0, q, q, 0.0, 0.0]
model['flowRate_time'] = [0.0,
                          t_init,
                          t_init+t_ramp,
                          t_init+t_ramp+t_pump,
                          t_init+(2*t_ramp)+t_pump,
                          1e9]

# Write the tables
write_table_files(model, output_path)


# Plot results
# plt.figure()
# plt.plot(model['sigma_xx']*-1e-6, z, label='sigma_xx')
# plt.plot(model['sigma_yy']*-1e-6, z, label='sigma_yy')
# plt.plot(model['sigma_zz']*-1e-6, z, label='sigma_zz')
# plt.plot(model['porePressure']*1e-6, z, label='P')
# plt.xlabel('Stress (MPa)')
# plt.ylabel('Depth (m)')
# plt.legend()

# plt.figure()
# plt.plot(model['nu'], z)
# plt.xlabel('Poissons Ratio')
# plt.ylabel('Depth (m)')

# plt.show()
