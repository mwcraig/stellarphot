from importlib import import_module

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.table import QTable
from astropy.utils.data import get_pkg_data_path

from stellarphot import SourceListData
from stellarphot.photometry import compute_fwhm, fast_fwhm_from_image, source_detection
from stellarphot.photometry.tests.fake_image import FakeCCDImage, FakeImage
from stellarphot.settings.models import FwhmMethods

# Make sure the tests are deterministic by using a random seed
SEED = 5432985


@pytest.mark.parametrize("units", [u.pixel, None])
def test_compute_fwhm(units):
    fake_image = FakeImage(seed=SEED)
    sources = fake_image.sources
    if units is not None:
        # It turns out having a unit on a column is not the same as
        # things in the column having units. The construct below ensures
        # that the source table values have units.
        # Do not try: sources['x_mean'] = sources['x_mean'] * units
        # Turns out individual values do NOT have units in that case.
        sources["x_mean"] = [v * units for v in sources["x_mean"]]
        sources["y_mean"] = [v * units for v in sources["y_mean"]]

    fwhm_x, fwhm_y = compute_fwhm(
        fake_image.image, sources, x_column="x_mean", y_column="y_mean"
    )

    expected_fwhm = np.array(sources["x_stddev"] * gaussian_sigma_to_fwhm)
    assert np.allclose(fwhm_x, expected_fwhm, rtol=1e-2)


def _image_with_bad_pixel(mask_by_nan):
    """
    Make a fake image in which the pixel at the center of the first source is
    bad, either because it is NaN or because it is masked.

    Parameters
    ----------
    mask_by_nan : bool
        If ``True``, mark the bad pixel by setting it to NaN in a plain numpy
        array. If ``False``, return a `~astropy.nddata.CCDData` whose mask is
        `True` at that pixel.

    Returns
    -------
    image : `numpy.ndarray` or `astropy.nddata.CCDData`
        The image with one bad pixel.
    sources : `astropy.table.Table`
        The table of sources in the image.
    bad_pixel : tuple of int
        The ``(x, y)`` position of the bad pixel.
    """
    fake_image = FakeImage(seed=SEED)
    sources = fake_image.sources
    x, y = sources["x_mean"].astype(int)[0], sources["y_mean"].astype(int)[0]
    image = fake_image.image.copy()

    if mask_by_nan:
        # Note the usual row/column swap when going to x/y coordinates.
        image[y, x] = np.nan
    else:
        image = CCDData(image, unit=u.adu, mask=np.zeros_like(image, dtype=bool))
        image.mask[y, x] = True

    return image, sources, (x, y)


def _expected_fwhm(sources):
    return np.array(sources["x_stddev"] * gaussian_sigma_to_fwhm)


@pytest.mark.parametrize("mask_by_nan", [True, False])
def test_compute_fwhm_with_missing_data(mask_by_nan):
    # Regression test for https://github.com/feder-observatory/stellarphot/issues/161
    # We should be able to find FWHM for a source even with NaNs in the image.
    image, sources, _ = _image_with_bad_pixel(mask_by_nan)

    fwhm_x, _ = compute_fwhm(
        image, sources, x_column="x_mean", y_column="y_mean", fit_method=FwhmMethods.FIT
    )

    assert np.allclose(fwhm_x, _expected_fwhm(sources), rtol=1e-2)


def test_compute_fwhm_does_not_zero_fill_nans(mocker):
    # Regression test for https://github.com/feder-observatory/stellarphot/issues/641
    # compute_fwhm used to replace NaN pixels with zero before calling
    # photutils' fit_fwhm, as a workaround for photutils issue #2029. That was
    # fixed in photutils 2.3, which is now the minimum version, so the NaNs
    # should reach fit_fwhm untouched (masked, but not overwritten).
    image, sources, (x, y) = _image_with_bad_pixel(mask_by_nan=True)

    # The function under test uses fit_fwhm from its own module namespace, and
    # that name is shadowed in the stellarphot.photometry namespace by the
    # source_detection function, so get at the module explicitly.
    detection_module = import_module("stellarphot.photometry.source_detection")
    spy = mocker.spy(detection_module, "fit_fwhm")

    fwhm_x, _ = compute_fwhm(
        image, sources, x_column="x_mean", y_column="y_mean", fit_method=FwhmMethods.FIT
    )

    # The first source is the one whose central pixel is NaN.
    first_call = spy.call_args_list[0]
    data_passed = first_call.args[0]
    x_cutout, y_cutout = first_call.kwargs["xypos"]

    # The NaN made it to photutils, and it is still at the position of the star.
    nan_positions = np.argwhere(np.isnan(data_passed))
    assert len(nan_positions) == 1
    nan_y, nan_x = nan_positions[0]
    assert abs(nan_x - x_cutout) <= 1
    assert abs(nan_y - y_cutout) <= 1

    # The NaN is masked in the call, and the FWHM is still measured correctly.
    assert first_call.kwargs["mask"][nan_y, nan_x]
    assert np.allclose(fwhm_x, _expected_fwhm(sources), rtol=1e-2)


@pytest.mark.parametrize("image_is_ccd", [True, False])
@pytest.mark.parametrize("find_fwhm", [True, False])
@pytest.mark.parametrize("stddev_input", [True, False])
@pytest.mark.parametrize("spp_unit", [None, u.adu / u.pixel])
@pytest.mark.parametrize("provide_spp_input", [True, False, None])
def test_detect_source_number_location(
    provide_spp_input, spp_unit, stddev_input, find_fwhm, image_is_ccd
):
    """
    Make sure we detect the sources in the input table....
    """
    # Skip some cmobinations we really don't need to test
    if not provide_spp_input and spp_unit is not None:
        pytest.skip("No point in testing this combination")
    if image_is_ccd:
        fake_image = FakeCCDImage(seed=SEED)
        image = fake_image
    else:
        fake_image = FakeImage(seed=SEED)
        image = fake_image.image

    sources = QTable(
        fake_image.sources,
        units={
            "x_mean": u.pixel,
            "y_mean": u.pixel,
            "x_stddev": u.pixel,
            "y_stddev": u.pixel,
        },
    )
    # print(sources)
    if provide_spp_input:
        # Pass only one value for the sky background for source detection
        sky_per_pix = sources["sky_per_pix_avg"].mean()
        sky_per_pix *= spp_unit if spp_unit else 1
    else:
        sky_per_pix = None
        print(sources["sky_per_pix_avg"].mean())

    stddev = fake_image.noise_dev if stddev_input else None

    fwhm = 2 * sources["x_stddev"].mean()

    found_sources = source_detection(
        image,
        fwhm=fwhm,
        find_fwhm=find_fwhm,
        threshold=10,
        sky_per_pix_avg=sky_per_pix,
        stddev=stddev,
        verbose=True,  # Run verbose=True on this test to cover all verbose lines
    )
    # Sort by flux so we can reliably match them
    sources.sort("amplitude")
    found_sources.sort("flux")

    # Do we have the right number of sources?
    assert len(sources) == len(found_sources)

    for inp, out in zip(sources, found_sources, strict=True):
        # Do the positions match?
        np.testing.assert_allclose(out["xcenter"], inp["x_mean"], rtol=1e-5, atol=0.05)
        np.testing.assert_allclose(out["ycenter"], inp["y_mean"], rtol=1e-5, atol=0.05)
        if find_fwhm:
            np.testing.assert_allclose(
                gaussian_sigma_to_fwhm * (inp["x_stddev"] + inp["y_stddev"]) / 2,
                out["width"],
                rtol=1e-5,
                atol=0.05,
            )


def test_detect_source_with_padding():
    """
    Make sure we detect the sources in the input table....
    """
    fake_image = FakeImage(seed=SEED)
    sources = QTable(
        fake_image.sources,
        units={
            "x_mean": u.pixel,
            "y_mean": u.pixel,
            "x_stddev": u.pixel,
            "y_stddev": u.pixel,
        },
    )
    # Pass only one value for the sky background for source detection
    sky_per_pix = sources["sky_per_pix_avg"].mean()
    # Padding was chosen to be large enough to ensure that one of the sources in
    # test_sources.csv would land too close to the edge of the image.
    found_sources = source_detection(
        fake_image.image,
        fwhm=2 * sources["x_stddev"].mean(),
        threshold=10,
        sky_per_pix_avg=sky_per_pix,
        padding=95,
    )

    # Did we drop one source because it was too close to the edge?
    assert len(sources) - 1 == len(found_sources)


def test_detect_source_bad_input():
    with pytest.raises(ValueError, match="ccd must be a numpy array or CCDData object"):
        source_detection(None)


@pytest.mark.parametrize(
    "fit_method", [FwhmMethods.FIT, FwhmMethods.PROFILE, FwhmMethods.MOMENTS]
)
def test_fwhm_computation(fit_method):
    # Regression test for #490, in which FWHM computation is incorrect
    # because the image hasn't been background subtracted.
    ccd_file = get_pkg_data_path("data/cutout_for_fwhm_test.fits")
    source_list_file = get_pkg_data_path("data/source_list_for_fwhm_test.ecsv")

    ccd = CCDData.read(ccd_file, unit=u.adu)
    source_list = SourceListData.read(source_list_file)
    # Value below is from the image the cutout was taken from
    source_list["sky_per_pix_avg"] = np.array([83.69])

    fwhm_x, fwhm_y = compute_fwhm(
        ccd,
        source_list,
        fwhm_estimate=6.85179,
        x_column="xcenter",
        y_column="ycenter",
        sky_per_pix_column="sky_per_pix_avg",
        fit_method=fit_method,
    )

    avg_fwhm = np.mean([fwhm_x, fwhm_y])
    if fit_method == FwhmMethods.MOMENTS:
        assert avg_fwhm > 7
    else:
        assert np.isclose(avg_fwhm, 6.6, rtol=0.1)


def test_compute_fwhm_input_options():
    with pytest.raises(ValueError, match="Cannot specify both "):
        compute_fwhm(
            None, None, sky_per_pix_avg=10, sky_per_pix_column="sky_per_pix_avg"
        )

    # Make a table with a single column
    table = QTable({"sky_per_pix_avg": [10]})
    spp_column = "foo"
    with pytest.raises(
        ValueError, match=f"Column {spp_column} not found in sources table"
    ):
        compute_fwhm(None, table, sky_per_pix_column=spp_column)

    fake_image = FakeImage(seed=SEED)
    sources = fake_image.sources

    fit_method = "foo"
    with pytest.raises(ValueError, match=f"Unknown fit method: {fit_method}"):
        compute_fwhm(
            fake_image.image,
            sources,
            fit_method=fit_method,
            x_column="x_mean",
            y_column="y_mean",
        )


@pytest.mark.parametrize(
    "aggregate_by,n_bright_expect",
    (
        [None, 20],
        ["mean", 1],
        ["median", 1],
    ),
)
def test_fast_fwhm_from_image(aggregate_by, n_bright_expect):
    # This should maybe be split into a few tests. Putting that aside, this test
    # checks that
    # + the FWHM estimate is correct
    # + changing the max_adu changes the sources that are used
    expected_fwhm = 5.5
    fake_image = FakeImage(seed=SEED, fwhm=expected_fwhm, n_repeats_per_side=5)
    n_bright = 20
    fwhms = fast_fwhm_from_image(
        fake_image.image,
        5,
        noise=fake_image.noise_dev,
        n_brightest_sources=n_bright,
        max_adu=65000,
        aggregate_by=aggregate_by,
    )

    # Make sure we got the right number of sources
    assert len(np.atleast_1d(fwhms)) == n_bright_expect

    if aggregate_by is None:
        # Check that the FWHM values are reasonable
        assert np.allclose(fwhms, fake_image.input_fwhm, rtol=0.01)
    else:
        assert fwhms == pytest.approx(fake_image.input_fwhm, rel=0.01)

    # Make sure some of the stars are over the max_adu and check that
    # we get a different list of sources.
    max_adu = fake_image.sources["amplitude"].max() * 0.6
    fwhms_some_max_out = fast_fwhm_from_image(
        fake_image.image,
        5,
        noise=fake_image.noise_dev,
        n_brightest_sources=n_bright,
        max_adu=max_adu,
        aggregate_by=aggregate_by,
    )

    if aggregate_by is None:
        # Make sure we got different sources
        assert fwhms.mean() != fwhms_some_max_out.mean()
    else:
        assert fwhms != fwhms_some_max_out


def test_fast_fwhm_from_image_bad_aggregate():
    # check that providing a bad aggregate_by method raises an error
    expected_fwhm = 5.5
    fake_image = FakeImage(seed=SEED, fwhm=expected_fwhm, n_repeats_per_side=5)
    n_bright = 20
    with pytest.raises(ValueError, match="Unknown aggregate_by method"):
        fast_fwhm_from_image(
            fake_image.image,
            5,
            noise=fake_image.noise_dev,
            n_brightest_sources=n_bright,
            max_adu=65000,
            aggregate_by="foo",
        )


def test_fast_fwhm_from_image_does_not_mutate_input_mask():
    # Regression test for a latent bug in fast_fwhm_from_image: the function
    # used to do `mask |= data > max_adu` on the mask pulled directly off of
    # the input CCDData, which mutates the caller's mask array in place.
    expected_fwhm = 5.5
    fake_image = FakeCCDImage(seed=SEED, fwhm=expected_fwhm)

    # Build a mask with a mix of True and False entries that is unrelated to
    # which pixels are above max_adu.
    mask = np.zeros(fake_image.data.shape, dtype=bool)
    mask[0, 0] = True
    mask[0, 1] = True
    mask[10, 10] = True
    fake_image.mask = mask

    original_mask = mask.copy()

    # Set max_adu comfortably above every source's peak so source detection
    # and FWHM fitting behave normally, then poke in a single hot pixel,
    # away from any source and not already masked, that exceeds max_adu.
    # This guarantees the `mask |= data > max_adu` line actually finds a
    # new pixel to mask, regardless of the random noise realization.
    max_adu = 5000.0
    assert fake_image.data.max() < max_adu
    fake_image.data[5, 5] = max_adu + 1000
    assert not fake_image.mask[5, 5]

    fast_fwhm_from_image(
        fake_image,
        5,
        noise=fake_image.noise_dev,
        max_adu=max_adu,
        aggregate_by=None,
    )

    np.testing.assert_array_equal(fake_image.mask, original_mask)


def test_block_center_to_pixel():
    # Regression test for the block-center off-by-0.5 bug in
    # fast_fwhm_from_image (stellarphot issue #606). A block of size
    # ``block_size`` at reduced-image index ``block_index`` covers original
    # image pixels ``block_index * block_size`` through
    # ``(block_index + 1) * block_size - 1`` inclusive, so the center of
    # that block, in original-image pixel coordinates, is
    # ``block_size * (block_index + 0.5) - 0.5``, NOT
    # ``block_size * (block_index + 0.5)``.
    from stellarphot.photometry.source_detection import _block_center_to_pixel

    for block_size in [1, 4, 8]:
        for block_index in [0, 1, 5, 12]:
            original_pixels = np.arange(
                block_index * block_size, (block_index + 1) * block_size
            )
            expected_center = original_pixels.mean()
            assert _block_center_to_pixel(block_index, block_size) == pytest.approx(
                expected_center
            )
