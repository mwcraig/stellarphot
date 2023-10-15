from astropy.table import Table, join
import numpy as np


def broeg_weights2(mags,
                   errors,
                   max_iterations=10,
                   convergence_threshold=0.0001,
                   error_factor=1.0,
                   sigma_nI=0.0,):
    """
    Calculate the Broeg weights for a given set of instrumental magnitudes and errors.

    Parameters
    ----------

    mags : array-like
        The instrumental magnitudes.

    errors : array-like
        The instrumental magnitude errors.

    max_iterations : int, optional
        The maximum number of iterations to perform.

    convergence_threshold : float, optional
        The threshold of difference in weights at which to stop iterating.

    error_factor : float, optional
        The factor by which to multiply the errors to get the initial weights.

    sigma_nI : float, optional
        The additive term to add to the errors to get the initial weights.

    Returns
    -------
    array-like
        The Broeg weights.

    Notes
    -----

    This is an implementation of the algorithm described in C. Broeg at al's
    algorithm ('A new algorithm for differential photometry:
    computing an optimum artificial comparison
    star', 2005, http://adsabs.harvard.edu/abs/2005AN....326..134B).
    """
    weight_vec = 1.0 / (error_factor * errors.mean(axis=1)**2 + sigma_nI**2)
    n_stars = mags.shape[0]

    weight_vecs = [weight_vec / weight_vec.sum()]
    light_curves = []
    for i in range(max_iterations):
        if i > 2:
            # Check for convergence
            if np.abs(weight_vecs[-1] - weight_vecs[-2]).max() < convergence_threshold:
                weight_vecs = weight_vecs[:-1]
                break

        # Make weights into an array
        weights = np.array(
            [
                weight_vec for _ in range(n_stars)
            ]
        )

        # Diagonal entries are zero because a star cannot be compared
        # to itself.
        weights[np.diag_indices(n_stars)] = 0.0

        # Normalize the weights -- note the normalization factor
        # is different for each star.
        weights_norm_factor = 1 / weights.sum(axis=1)

        # Create a diagonal matrix of the normalization factors
        weights_norm_matrix = np.zeros([n_stars, n_stars])
        weights_norm_matrix[np.diag_indices(n_stars)] = weights_norm_factor

        # With the previous definitions, the comparison magnitudes, i.e.
        # the magnitudes of the artificial comparison star, are given by
        # the following matrix multiplication:
        comp_mags = weights_norm_matrix @ weights @ mags

        # Find the differential magnitude for each star...
        delta_mags = mags - comp_mags
        # ...and save it.
        light_curves.append(delta_mags)

        # Calculate the new weights, which are the inverse of the
        # variance of the differential magnitudes.
        weight_vec = 1 / delta_mags.std(axis=1)**2
        weight_vecs.append(weight_vec / weight_vec.sum())

    return weight_vecs, light_curves
