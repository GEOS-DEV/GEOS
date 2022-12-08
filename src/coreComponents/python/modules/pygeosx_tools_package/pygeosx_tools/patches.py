

def patch_numpy_where():
    """
    Some packages (numpy, etc.) use np.where in a way that isn't
    compliant with the API.  In most cases this is fine, but when the
    geosx environment is initialized this can cause critical errors with
    nan inputs.  This patch checks and updates the inputs to where
    """
    import numpy as np
    np.where_original_fn = np.where
    print('patching np.where', flush=True)


    def variable_as_array_type(x, target_shape):
        if x is None:
            return

        if isinstance(x, (np.ndarray, list)):
            return x

        y = np.empty(target_shape)
        y[...] = x
        return y

    def flexible_where(condition, x=None, y=None):
        s = np.shape(condition)
        x = variable_as_array_type(x, s)
        y = variable_as_array_type(y, s)
        return np.where_original_fn(condition, x, y)

    np.where = flexible_where



patch_numpy_where()
