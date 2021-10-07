
import numpy as np


def write_table_files(tables, output_path='./', output_format='%1.5e'):
  for k in tables.keys():
    tmp = np.reshape(tables[k], (-1), order='F')
    np.savetxt('%s/%s.csv' % (output_path, k), tmp, fmt=output_format, delimiter=',')


# Config
zo = 500.0
rho = 2650.0
stress_std = 0.5e6
modulus_std = 0.5e9
shear_modulus_mean = 8.3e+09
bulk_modulus_mean = 11.9e+10

# x = [0.0, 5.0, 10.0]
# y = [0.0, 5.0, 10.0]
# z = [0.0, 5.0, 10.0]
x = [0.0]
y = [0.0]
z = np.linspace(0, 10, 4)
np.random.seed(12345)

# Build table axes
NZ = len(z)
tables = {'x': x, 'y': y, 'z': z}
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Generate stress model with vertical variations
tables['sigma_zz'] = -9.81 * rho * (Z + zo)
tables['sigma_xx'] = 0.5 * tables['sigma_zz'] + 2e6
tables['sigma_yy'] = 0.5 * tables['sigma_zz'] + 3e6
# for k in ['sigma_xx', 'sigma_yy']:
#   for ii in range(0, NZ):
#     tables[k][:, :, ii] += stress_std * np.random.randn()

# Generate random shear, bulk modulus
tables['bulkModulus'] = bulk_modulus_mean + (modulus_std * np.random.randn(*np.shape(X)))
tables['shearModulus'] = shear_modulus_mean + (modulus_std * np.random.randn(*np.shape(X)))

# Write tables
write_table_files(tables)

