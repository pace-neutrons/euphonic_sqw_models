import os

import pytest
import numpy as np
import numpy.testing as npt
import scipy.io
import euphonic
from euphonic import ureg

import euphonic_sqw_models

def get_abspath(filename, sub_dir):
    test_dir = os.path.dirname(os.path.realpath(__file__))
    return test_dir + os.path.sep + sub_dir + os.path.sep + filename

temp = [300]
quartz = ['quartz', euphonic.ForceConstants.from_castep,
          {'filename': get_abspath('quartz.castep_bin', 'input')}]
nacl =  ['nacl', euphonic.ForceConstants.from_phonopy,
         {'path': get_abspath('NaCl', 'input')}]
nacl_json =  ['nacl', euphonic.ForceConstants.from_json_file,
              {'filename': get_abspath('nacl_force_constants.json', 'input')}]
dw_grid = [None, [6,6,6]]
bose = [None, False]
negative_e = [None, True]
conversion_mat = [None, (1./2)*np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])]
lim = [None, 1e2]  # Units: mbarn


qpts = [[0.0,  0.0,  0.0],
        [0.1,  0.2,  0.3],
        [0.4,  0.5,  0.0],
        [0.6,  0.0,  0.7],
        [0.0,  0.8,  0.9],
       [-0.5,  0.0,  0.0],
        [0.0, -0.5,  0.0],
        [0.0,  0.0, -0.5],
        [1.0, -1.0, -1.0]]
qpts = np.array(qpts)
scattering_lengths = {'La': 8.24, 'Zr': 7.16, 'O': 5.803,
                      'Si': 4.1491, 'Na': 3.63, 'Cl': 9.577}
scattering_lengths_quantity = {'La': 8.24*ureg('fm'), 'Zr': 7.16*ureg('fm'),
                               'O': 5.803*ureg('fm'), 'Si': 4.1491*ureg('fm'),
                               'Na': 3.63*ureg('fm'), 'Cl': 9.577*ureg('fm')}
iscale = 1.0
freq_scale = 1.0
pars = [iscale, freq_scale]


def get_test_opts():
    opts = []
    for tt in temp:
        for dw in dw_grid:
            for bos in bose:
                for ne in negative_e:
                    for cmat in conversion_mat:
                        for lm in lim:
                            opts.append({
                                'temperature': tt,
                                'debye_waller_grid': dw,
                                'bose': bos,
                                'negative_e': ne,
                                'conversion_mat': cmat,
                                'lim': lm})
    return opts


def get_expected_output_dir():
    return 'expected_output'


def get_expected_output_filename(material_name, opts):
    fname = f'{material_name}'
    if opts.get('temperature', None) is not None:
        if hasattr(opts.get('temperature'), 'units'):
            fname += f'_T{opts["temperature"].magnitude}'
        else:
            fname += f'_T{opts["temperature"]}'
    if opts.get('debye_waller_grid', None) is not None:
        fname += '_dw' + ''.join([str(v) for v in opts['debye_waller_grid']])
    if opts.get('bose', None) is not None:
        fname += '_bose' + str(opts['bose']).lower()
    if opts.get('negative_e', None) is not None:
        fname += '_negative_e' + str(opts['negative_e']).lower()
    if opts.get('conversion_mat') is not None:
        fname += '_conv_mat_det' + str(np.linalg.det(opts['conversion_mat']))
    if opts.get('lim', None) is not None:
        # Account for different str conversion of positive exponent
        # sci notation - e.g. in Matlab string(1e2) = "100" but in
        # Python str(1e2) = "100.0" unless explicitly cast to int
        if opts['lim']%1 == 0:  # fractional part is zero, is int
            strlim = str(int(opts['lim']))
        else:
            strlim = str(opts['lim'])
        fname += '_lim' + strlim
    fname = fname.replace('.', 'p').replace('-', 'm')
    fname += '.mat'
    return get_abspath(fname, get_expected_output_dir());


def sum_degenerate_modes(w, sf):
    tol = 0.1;
    rows, cols = np.shape(w);
    summed_sf = np.zeros(np.shape(sf));
    for i in range(rows):
        sum_at = np.where(np.diff(w[i,:]) > tol)[0];
        x = np.zeros(cols);
        x[sum_at + 1] = 1;
        degenerate_modes = np.asarray(np.cumsum(x), dtype=int);
        summed = np.bincount(degenerate_modes, weights=sf[i,:]);
        summed_sf[i, :len(summed)] = summed;
    return summed_sf


def calculate_w_sf(material_opts, material_constructor, opt_dict,
                   sum_sf=True, hdisp_args=(), hdisp_kwargs={}):
    fc = material_constructor(**material_opts)
    opt_dict['asr'] = 'reciprocal'
    coherent_sqw = euphonic_sqw_models.CoherentCrystal(fc, **opt_dict)
    w, sf = coherent_sqw.horace_disp(
        qpts[:,0], qpts[:,1], qpts[:,2], *hdisp_args, **hdisp_kwargs)
    w = np.array(w).T
    sf = np.array(sf).T
    # Ignore gamma-point frequencies by setting to zero - their
    # values are unstable
    gamma_points = (np.sum(np.absolute(qpts - np.rint(qpts)), axis=-1)
                    < 1e-10)
    w[gamma_points, :3] = 0.
    # Ignore acoustic structure factors by setting to zero - their
    # values can be unstable at small frequencies
    sf[:,:3] = 0.
    if 'negative_e' in opt_dict and opt_dict['negative_e']:
        n = int(np.shape(sf)[1] / 2)
        sf[:,n:(n+3)] = 0
    if sum_sf:
        sf = sum_degenerate_modes(w, sf)
    return w, sf


def get_expected_w_sf(fname, sum_sf=True):
    ref_dat = scipy.io.loadmat(fname)
    # Was written from Matlab - convert first
    if ref_dat['expected_w'].dtype == object:
        expected_w = np.concatenate(ref_dat['expected_w'][0], axis=1)
        expected_sf = np.concatenate(ref_dat['expected_sf'][0], axis=1)
    else:
    # Was written from Python - use as is
        expected_w = ref_dat['expected_w']
        expected_sf = ref_dat['expected_sf']
    # Set gamma-point acoustic frequencies to zero
    gamma_points = (np.sum(np.absolute(qpts - np.rint(qpts)), axis=-1)
                    < 1e-10)
    expected_w[gamma_points, :3] = 0.
    # Set acoustic structure factors to zero
    expected_sf[:,:3] = 0
    if 'negative_e' in fname:
        n = int(np.shape(expected_sf)[1] / 2)
        expected_sf[:,n:(n+3)] = 0
    if sum_sf:
        expected_sf = sum_degenerate_modes(expected_w, expected_sf)
    return expected_w, expected_sf


def run_and_test_horace_disp(material, opt_dict, filename=None):
    material_name, material_constructor, material_opts = material
    if filename is None:
        filename = get_expected_output_filename(material_name, opt_dict)
    expected_w, expected_sf = get_expected_w_sf(filename)
    w, sf = calculate_w_sf(material_opts, material_constructor, opt_dict)
    npt.assert_allclose(w, expected_w, atol=1e-2, rtol=1e-5)
    npt.assert_allclose(sf, expected_sf, rtol=1e-2, atol=1e-2)


@pytest.mark.parametrize('material', [quartz])
@pytest.mark.parametrize('opt_dict', get_test_opts())
# The following options shouldn't change the result
@pytest.mark.parametrize('run_opts', [
    {'use_c': False, 'n_threads': 1, 'chunk': 5, 'dipole_parameter': 0.75,
     'scattering_lengths': scattering_lengths},
    {'use_c': True, 'n_threads': 1, 'chunk': 0,
     'scattering_lengths': scattering_lengths_quantity},
    {'use_c': True, 'n_threads': 2}])
def test_euphonic_sqw_models(material, opt_dict, run_opts):
    opt_dict.update(run_opts)
    # Don't pass None options
    opt_dict = {k: v for k, v in opt_dict.items() if v is not None}
    run_and_test_horace_disp(material, opt_dict)


@pytest.mark.phonopy_reader
@pytest.mark.parametrize('material', [nacl])
@pytest.mark.parametrize('opt_dict', get_test_opts())
def test_euphonic_sqw_models_phonopy(material, opt_dict):
    # Don't pass None options
    opt_dict = {k: v for k, v in opt_dict.items() if v is not None}
    run_and_test_horace_disp(material, opt_dict)


@pytest.mark.parametrize('material', [quartz])
def test_euphonic_sqw_models_defaults(material):
    filename = get_abspath(f'{material[0]}_defaults.mat',
                           get_expected_output_dir())
    run_and_test_horace_disp(material, {}, filename=filename)


@pytest.mark.phonopy_reader
@pytest.mark.parametrize('material', [nacl])
def test_euphonic_sqw_models_defaults_phonopy(material):
    filename = get_abspath(f'{material[0]}_defaults.mat',
                           get_expected_output_dir())
    run_and_test_horace_disp(material, {}, filename=filename)


@pytest.mark.parametrize('material', [quartz])
@pytest.mark.parametrize('temperature', [None, 0*ureg('K')])
@pytest.mark.parametrize('opt_dict', [{
    'debye_waller_grid': [6, 6, 6],
    'negative_e': True,
    'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                       [1, -1, 1],
                                       [1, 1, -1]])}])
def test_euphonic_sqw_models_temperatures(
        material, opt_dict, temperature):
    opt_dict.update({'temperature': temperature})
    run_and_test_horace_disp(material, opt_dict)


@pytest.mark.phonopy_reader
@pytest.mark.parametrize('material', [nacl])
@pytest.mark.parametrize('temperature', [None, 0*ureg('K')])
@pytest.mark.parametrize('opt_dict', [{
    'debye_waller_grid': [6, 6, 6],
    'negative_e': True,
    'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                       [1, -1, 1],
                                       [1, 1, -1]])}])
def test_euphonic_sqw_models_temperatures_phonopy_reader(
        material, opt_dict, temperature):
    opt_dict.update({'temperature': temperature})
    run_and_test_horace_disp(material, opt_dict)


@pytest.mark.parametrize('material', [quartz, nacl_json])
@pytest.mark.parametrize('opt_dict', [{
    'temperature': 300,
    'debye_waller_grid': [6, 6, 6],
    'negative_e': True,
    'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                       [1, -1, 1],
                                       [1, 1, -1]])}])
# Test a combination of keyword and positional arguments
# Is required for use with pace-python (Toby/Multifit
# allow positional only for fitting parameters)
@pytest.mark.parametrize('iscale, freqscale, args, kwargs',
        [(1.0, 1.3, (), {'frequency_scale': 1.3}),
         (1e-4, 1.0, (1e-4,), {}),
         (1.5e3, 0.4, (1.5e3, 0.4), {}),
         (2e3, 0.5, (), {'frequency_scale': 0.5, 'intensity_scale': 2e3}),
         (50, 0.9, (50,), {'frequency_scale': 0.9})])
def test_euphonic_sqw_models_pars(
        material, opt_dict, iscale, freqscale, args, kwargs):
    material_name, material_constructor, material_opts = material
    expected_w, expected_sf = get_expected_w_sf(
        get_expected_output_filename(
            material_name, opt_dict), sum_sf=False)
    w, sf = calculate_w_sf(material_opts, material_constructor, opt_dict,
                           sum_sf=False, hdisp_args=args, hdisp_kwargs=kwargs)
    # Sum sfs using the same frequencies - if expected_w and the scaled
    # w are used, this could result in expected_sf and sf
    # being summed differently
    sf_summed = sum_degenerate_modes(expected_w, sf)
    expected_sf_summed = sum_degenerate_modes(expected_w, expected_sf)
    # Check that pars have scaled w and sf as expected
    npt.assert_allclose(w, expected_w*freqscale, atol=1e-2*freqscale,
                        rtol=1e-5)
    npt.assert_allclose(sf_summed, iscale*expected_sf_summed,
                        rtol=1e-2, atol=1e-2*iscale)


@pytest.mark.parametrize('opt_dict', [{
    'temperature': 300*ureg('K'),
    'debye_waller_grid': [6, 6, 6],
    'negative_e': True,
    'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                       [1, -1, 1],
                                       [1, 1, -1]])}])
def test_old_behaviour_single_parameter_sets_intensity_scale(opt_dict):
    iscale = 1.5
    material_name, material_constructor, material_opts = quartz
    fc = material_constructor(**material_opts)
    opt_dict['asr'] = 'reciprocal'
    opt_dict['scattering_lengths'] = scattering_lengths
    coherent_sqw = euphonic_sqw_models.CoherentCrystal(fc, **opt_dict)
    # Pass a single positional argument for iscale - this is the old
    # syntax. Test that it hasn't broken, if it has we may need to raise
    # a deprecation warning
    w, sf = coherent_sqw.horace_disp(
        qpts[:,0], qpts[:,1], qpts[:,2], [iscale])
    w = np.array(w).T
    sf = np.array(sf).T
    # Ignore gamma-point frequencies by setting to zero - their
    # values are unstable
    gamma_points = (np.sum(np.absolute(qpts - np.rint(qpts)), axis=-1)
                    < 1e-10)
    w[gamma_points, :3] = 0.
    # Ignore acoustic structure factors by setting to zero - their
    # values can be unstable at small frequencies
    n = int(np.shape(sf)[1] / 2)
    sf[:,:3] = 0.
    sf[:,n:(n+3)] = 0

    expected_w, expected_sf = get_expected_w_sf(
        get_expected_output_filename(
            material_name, opt_dict), sum_sf=False)
    # Sum sfs using the same frequencies - if expected_w and the scaled
    # w are used, this could result in expected_sf and sf
    # being summed differently
    sf_summed = sum_degenerate_modes(expected_w, sf)
    expected_sf_summed = sum_degenerate_modes(expected_w, expected_sf)
    # Check that pars have scaled w and sf as expected
    npt.assert_allclose(w, expected_w, atol=1e-2, rtol=1e-5)
    npt.assert_allclose(sf_summed, iscale*expected_sf_summed,
                        rtol=1e-2, atol=1e-2)


# Units of sf have changed, test the calculated and old expected
# values are the same apart from a scale factor
@pytest.mark.parametrize('material', [quartz, nacl_json])
@pytest.mark.parametrize('opt_dict', [{
    'temperature': 300,
    'debye_waller_grid': [6, 6, 6],
    'negative_e': True,
    'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                       [1, -1, 1],
                                       [1, 1, -1]])}])
def test_sf_unit_change(material, opt_dict):
    material_name, material_constructor, material_opts = material
    fname = get_expected_output_filename(
            material_name, opt_dict)
    expected_w, expected_sf = get_expected_w_sf(os.path.join(
        os.path.dirname(fname),
        f'old_units_{os.path.basename(fname)}'))

    w, sf = calculate_w_sf(material_opts, material_constructor, opt_dict)

    npt.assert_allclose(w, expected_w, atol=1e-2, rtol=1e-5)

    # Change in units angstrom**2 -> mbarn = 1e11. Change from sf per
    # unit cell to sf per atom = 1/n_atoms. Change from relative to
    # absolute units means 1/2 factor is now included
    n_atoms = w.shape[1]/6 # n_atoms is n_modes/6 because of 'negative_e'
    unit_scale = 1e11/(2*n_atoms)
    npt.assert_allclose(sf, expected_sf*unit_scale, rtol=1e-2, atol=1e-2)

def test_invalid_fcs_raises_value_error():
    with pytest.raises(ValueError):
        euphonic_sqw_models.CoherentCrystal(None)

@pytest.mark.parametrize('kwarg', [{'debye_waller': np.random.rand(5, 3, 3)},
                                   {'debye_waller_grid': [2, 3]},
                                   {'conversion_mat': np.random.rand(2, 3)},
                                   {'temperature': -5*ureg('K')}])
def test_invalid_kwargs_raises_value_error(kwarg):
    fc = euphonic.ForceConstants.from_castep(
        get_abspath('quartz.castep_bin', 'input'))
    with pytest.raises(ValueError):
        euphonic_sqw_models.CoherentCrystal(fc, **kwarg)

