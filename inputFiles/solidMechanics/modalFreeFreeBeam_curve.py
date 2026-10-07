import numpy as np


def _shape(x, beta_length, length):
    """Euler-Bernoulli free-free mode shape, normalized so that its mean square over the beam is one."""
    bl = beta_length
    sigma = (np.cosh(bl) - np.cos(bl)) / (np.sinh(bl) - np.sin(bl))
    b = bl / length
    return np.cosh(b * x) + np.cos(b * x) - sigma * (np.sinh(b * x) + np.sin(b * x))


def _mode_shape(kwargs, field, beta_length):
    """Mass-normalized free-free bending mode of the beam of modalFreeFreeBeam_base.xml along its axis.

    The modes of a beam with a square cross-section come in pairs with the same frequency, so the polarization
    of the computed mode in the (y, z) plane is arbitrary, as is its sign. The polarization and the sign are
    taken from the computed mode: its component along the analytical shape. The check is on the shape and on
    the amplitude, which is set by the normalization phi^T M phi = 1, that is rho A integral(phi^2) = 1.
    """
    density = 7850.0
    length = 10.0
    side = 0.2
    area = side * side

    x = np.squeeze(kwargs[f'{field} ReferencePosition axis'][0, :, 0])
    computed = kwargs[f'{field} axis']
    shape = _shape(x, beta_length, length) / np.sqrt(density * area * length)

    # Polarization and sign: projection of the computed transverse displacement on the analytical shape
    direction = np.array([np.sum(computed[0, :, 1] * shape), np.sum(computed[0, :, 2] * shape)])
    direction /= np.linalg.norm(direction)

    expected = np.zeros(np.shape(computed))
    expected[:, :, 1] = direction[0] * shape
    expected[:, :, 2] = direction[1] * shape
    return expected


def first_bending(**kwargs):
    return _mode_shape(kwargs, 'modeShape7', 4.730040745)


def second_bending(**kwargs):
    return _mode_shape(kwargs, 'modeShape9', 7.853204624)
