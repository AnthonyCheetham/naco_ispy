import configparser
import os

# import dill as pickle
import numpy as np
from astropy import units as u
from astropy.io import fits

import graphic_contrast_lib

from dicpm.instrument_parameters import Instrument
from dicpm.embed_shell import ipsh
from dicpm.parameters import Reduction_parameters
from dicpm.reduction_wrapper import run_complete_reduction

""" Parameter setup for test """
used_instrument = Instrument(
    name='NACO',
    pixel_scale=u.pixel_scale(0.02719 * u.arcsec / u.pixel),
    telescope_diameter=8.2 * u.m)
# used_instrument = Instrument(
#     name='IRDIS',
#     pixel_scale=u.pixel_scale(0.01225 * u.arcsec / u.pixel),
#     telescope_diameter=8.2 * u.m)

reduction_parameters = Reduction_parameters(
    search_region=None,  # if None the search region is defined by the following two parameters, otherwise a boolean mask defining the search region can be given here
    search_region_inner_bound=2,  # inner radius in pixels
    search_region_outer_bound=40, # outer radius in pixels
    oversampling=1,
    data_auto_crop = True,
    data_crop_size=200,
    right_handed=False,  # Rotation durection (naco==False, sphere==true, lbt==??)
    # known_companion_position=None,  # [-33.1, -14.46] coordinates in pixels, to exclude known object from the model
    use_multiprocess=True,
    ncpus=6,
    prefix='',  # Add to file name
    result_folder='./DICPM/',
    # Ignore this part for now: for further developments
    reduce_single_position=False,
    inject_fake=False,
    true_position=None,  # np.array([20, 20]),
    guess_position=None,  # np.array([20, 20]),
    true_contrast=None,  # 5.981e-06,
    guess_contrast=None,  # 5.981e-06,
    number_of_pca_regressors=50,  # Not relevant for search, only for single reduction
    #####
    yx_anamorphism=[1., 1.],
    scaling='temp-median',  # Scaling of pca training set: alternative temp-quartile
    annulus_width=9,  # in pixels
    annulus_offset=0,
    autosize_psf_mask_in_lambda_over_d=True,  # radius of the unsaturated PSF I am using as a model, in pixels
    psf_mask_size_in_lambda_over_d=2.4,
    signal_mask_psf_size=21,  # if autosize None, fix size in pixels, must be >= reduction mask size
    reduction_mask_psf_size=21,  # if autosize None, fix size in pixels
    add_radial_regressors=True,  # Add regressors at same azimuth as search location, but displace by following values
    radial_separation_from_source=[-8, 6],  # -8, 6  -5, 3
    include_opposite_regressors=True,  # Add region opposite reduction region to regressor pool
    # This could be interesting to have a look at (but rememebr to switch off 'use_multiprocess' if using this)
    plot_all_diagnostics=False,
    verbose=False)

# Number of components in fractions of the maximum number, if more than one number is given
# The algorithm will loop over them
number_of_components_fraction = [0.3]
contrast_curve = True
# READ DATA
# I need a datacube
# Datacube, shape=[optional dimension for wavelenghts][#number of frames][data:data]
data_full = fits.getdata('master_cube_PCA.fits')
pa = np.loadtxt('parallactic_angle.txt')  # parallactic angle
wavelengths = np.array([3.8]) * u.micron  # wavelength of data in microns
bad_frames = None  # array with indices of bad frames to be removed
# Shape for flux psf: shape=[optional dimension for wavelength][data:data] (no timeseries)
flux_psf_full = fits.getdata('flux.fits')  # unsaturated and scaled PSF


# Scale PSF
flux_hdr = fits.getheader('flux.fits')
cube_hdr = fits.getheader('master_cube_PCA.fits')
flux_fact = graphic_contrast_lib.correct_flux_frame_naco(cube_hdr,flux_hdr)
print('Scaling flux frame by '+str(flux_fact))
flux_psf_full *= flux_fact

# Make the PSF have an odd number of pixels
flux_psf_full = flux_psf_full[1:,1:]

# Average before and after (average psf if I have more images for PSF), normally not useful for me

# Boolean image for each wavelength with True=bad pixel and False=good pixel
# useful when running reduction on raw unshifted data
bad_pixel_mask_full = None

# Centering
xy_image_centers = None  # if running reduction on un-centered data, provide xy position of center (in uncropped frame)

# Waffle amplitudes
amplitude_modulation_full = None  # modulation of model compared to median flux
contrast_map_full = None

# Example for small circular region at planet position
# reduction_parameters.search_region = regressor_selection.make_signal_mask(
#     yx_dim=(reduction_parameters.data_crop_size, reduction_parameters.data_crop_size),
#     pos_yx=[-33.1, -14.46],
#     mask_radius=15,
#     oversampling=reduction_parameters.oversampling)


run_complete_reduction(
    data_full=data_full,
    flux_psf_full=flux_psf_full,
    pa=pa,
    bad_frames=bad_frames,
    wavelengths=wavelengths,
    instrument=used_instrument,
    number_of_components_fraction=number_of_components_fraction,
    reduction_parameters=reduction_parameters,
    bad_pixel_mask_full=bad_pixel_mask_full,
    xy_image_centers=xy_image_centers,
    amplitude_modulation_full=amplitude_modulation_full,
    contrast_map_full=contrast_map_full,
    contrast_curve=contrast_curve)
