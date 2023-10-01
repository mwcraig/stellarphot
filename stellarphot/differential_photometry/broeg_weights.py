from astropy.table import Table, join
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

    use_phot = photometry[photometry['filter'] == filter]

    # Get number of stars in the photometry table. Do this after selecting the desired filter
    # in case there are observations in filters other than the one the user wants to use.
    star_ids = sorted(set(use_phot['star_id']))
    n_stars = len(star_ids)

    if n_stars < 3:
        raise ValueError('Photometry table must contain at least three stars')

    # Estimate the initial error by aggregating the error for each star over all
    # of the measurements of the error.
    use_phot_grouped = use_phot.group_by('star_id')

    avg_phot_grouped = use_phot_grouped.groups.aggregate(np.nanmean)

    # import pdb; pdb.set_trace()

    # Calculate the initial Broeg weights from the instrumental errors
    weights = 1.0 / avg_phot_grouped['mag_error']**2

    # Weights need to be normalized
    weights /= np.sum(weights)

    weights_list = [weights]
    light_curves = []

    # This will be slow, but easy to implement for now
    for i in range(max_iterations):
        weights = np.zeros_like(weights)
        for idx, a_star in enumerate(star_ids):
            # Calculate the weighted mean magnitude of the comparison stars, excluding
            # the current star
            weight_table = Table([star_ids, weights_list[-1]], names=['star_id', 'weight'])
            comp_mag = calc_comp_mag(use_phot, weight_table, exclude_star_id=a_star)

            this_star = use_phot['star_id'] == a_star
            # Calculate the difference between the two
            delta_mag = use_phot[this_star]['mag_inst'] - comp_mag

            # Calculate the new weight for the current star
            if np.isnan(delta_mag).all() or np.nanstd(delta_mag) == 0:
                print(f"Setting weight to 0 for star {a_star}")
                weights[idx] = 0.0
            else:
                weights[idx] = 1.0 / np.nanstd(delta_mag)**2
        weights_list.append(weights / np.sum(weights))
    return star_ids, weights_list


def calc_comp_mag(photometry, weights_table, exclude_star_id=None):
    """
    Calculate the weighted mean magnitude of the comparison stars.

    Parameters
    ----------

    photometry : `astropy.table.Table`
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
        weights_table = weights_table.copy()
        # Set the weight of the excluded star to zero...
        excluded_star_index = weights_table['star_id'] == exclude_star_id
        weights_table['weight'][excluded_star_index] = 0.0
        # ...and renormalize the weights
        weights_table['weight'] /= np.sum(weights_table['weight'])

    joined = join(photometry, weights_table, keys='star_id', join_type='left')
    joined['weighted_mag'] = joined['mag_inst'] * joined['weight']
    joined = joined['BJD', 'weighted_mag']
    agged = joined.group_by('BJD').groups.aggregate(np.nansum)
    # Calculate the weighted mean magnitude
    return agged['weighted_mag']


def broeg_weights2(mags, errors, max_iterations=10, convergence_threshold=0.0001):
    """
    Calculate the Broeg weights for a given set of instrumental magnitudes and errors.

    Parameters
    ----------

    mags : array-like
        The instrumental magnitudes.

    errors : array-like
        The instrumental magnitude errors.

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
    heavily on the LEMON package at
    """
    weight_vec = 1.0 / (errors.mean(axis=1)**2)
    n_stars = mags.shape[0]

    weight_vecs = [weight_vec / weight_vec.sum()]
    light_curves = []
    for i in range(max_iterations):
        if i > 2:
            # Check for convergence
            if np.abs(weight_vecs[-1] - weight_vecs[-2]).max() < convergence_threshold:
                weight_vecs = weight_vecs[:-1]
                break

        weights = np.array(
            [
                weight_vec for _ in range(n_stars)
            ]
        )
        weights[np.diag_indices(n_stars)] = 0.0
        weights_norm_factor = 1 / weights.sum(axis=1)
        weights_norm_matrix = np.zeros([n_stars, n_stars])
        weights_norm_matrix[np.diag_indices(n_stars)] = weights_norm_factor
        comp_mags = weights_norm_matrix @ weights @ mags
        delta_mags = mags - comp_mags
        light_curves.append(delta_mags)
        # print(f"Iteration {i}: {weight_vec}")
        # print(f"\tNormalized weights: {weights_norm_matrix @ weights}")
        # print(f"Delta mags: {delta_mags}")

        weight_vec = 1 / delta_mags.std(axis=1)**2
        weight_vecs.append(weight_vec / weight_vec.sum())

    return weight_vecs, light_curves
