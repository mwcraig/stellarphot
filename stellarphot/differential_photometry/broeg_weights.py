import numpy as np


def broeg_weights(photometry, filter=None, max_iterations=100, convergence_threshold=0.001):
    """
    Calculate the Broeg weights for a given set of instrumental magnitudes and errors.

    Parameters
    ----------

    photometry : `astropy.table.Table`
        The photometry table.

    filter : str, optional
        The filter to use. If not specified, and only one filter is present, that
        filter is used. If more than one filter is present an error is raised.

    max_iterations : int, optional
        The maximum number of iterations to perform.

    convergence_threshold : float, optional
        The convergence threshold for the iterative algorithm.

    Returns
    -------
    array-like
        The Broeg weights.

    Notes
    -----

    This implementation of the algorithm described in C. Broeg's
    algorithm ('A new algorithm for differential photometry:
    computing an optimum artificial comparison
    star', 2005, http://adsabs.harvard.edu/abs/2005AN....326..134B) is based
    heavily on the LEMON package at https://github.com/vterron/lemon/tree/master

    Specifically, see line 474 of diffphot.py in that package.
    """
    # Rather than following the LEMON approach of calculating intial weights from
    # instrumental magnitudes, we use the instrumental magnitude errors directly,
    # as Broeg suggests in the paper.

    # Get number of stars in the photometry table
    star_ids = set(photometry['star_id'])
    n_stars = len(star_ids)

    # What filters are present in the data?
    filters = set(photometry['filter'])

    # If a filter is specified, check that it is present in the data
    if filter is not None:
        if filter not in filters:
            raise ValueError('Filter {} not present in photometry table'.format(filter))

    else:
        # If no filter is specified, and only one filter is present, use that filter
        if len(filters) == 1:
            filter = filters.pop()
        else:
            raise ValueError('More than one filter present in photometry table. You must specify a filter.')

    if n_stars < 3:
        raise ValueError('Photometry table must contain at least three stars')

    use_phot = photometry[photometry['filter'] == filter]
    # Calculate the initial Broeg weights from the instrumental errors
    weights = 1.0 / use_phot['mag_error']**2

    # Weights need to be normalized
    weights /= np.sum(weights)

    # This will be slow, but easy to implement for now


    return weights


def calc_comp_mag(mags, weights, exclude_star_id=None):
    """
    Calculate the weighted mean magnitude of the comparison stars.

    Parameters
    ----------

    mags : array-like
        The instrumental magnitudes of the comparison stars.

    weights : array-like
        The Broeg weights.

    exclude_star_id : str, optional
        The ID of a star to exclude from the calculation.

    Returns
    -------
    float
        The weighted mean magnitude of the comparison stars.
    """
    # If a star is to be excluded, set its weight to zero
    if exclude_star_id is not None:
        weights = weights.copy()
        # Set the weight of the excluded star to zero...
        weights[exclude_star_id] = 0.0
        # ...and renormalize the weights
        weights /= np.sum(weights)

    # Calculate the weighted mean magnitude
    return np.sum(mags * weights)
