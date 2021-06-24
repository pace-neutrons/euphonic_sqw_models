import os

import pytest
import numpy as np
import numpy.testing as npt
import scipy.io
import euphonic

import euphonic_sqw_models

def get_abspath(filename, sub_dir):
    test_dir = os.path.dirname(os.path.realpath(__file__))
    return test_dir + os.path.sep + sub_dir + os.path.sep + filename

temp = [300]
materials = [['quartz', euphonic.ForceConstants.from_castep,
              {'filename': get_abspath('quartz.castep_bin', 'input')}],
             ['nacl', euphonic.ForceConstants.from_phonopy,
              {'path': get_abspath('NaCl', 'input')}]]
dw_grid = [None, [6,6,6]]
bose = [None, False]
negative_e = [None, True]
conversion_mat = [None, (1./2)*np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])]
lim = [None, 1e2]

run_pars = [{'use_c': False, 'n_threads': 1, 'chunk': 5, 'dipole_parameter': 0.75},
            {'use_c': True, 'n_threads': 1, 'chunk': 0},
            {'use_c': True, 'n_threads': 2, 'chunk': 0}]

qpts = [[0.0,  0.0,  0.0],
        [0.1,  0.2,  0.3],
        [0.4,  0.5,  0.0],
        [0.6,  0.0,  0.7],
        [0.0,  0.8,  0.9],
       [-0.5,  0.0,  0.0],
        [0.0, -0.5,  0.0],
        [0.0,  0.0, -0.5],
        [1.0, -1.0, -1.0]]
scattering_lengths = {'La': 8.24, 'Zr': 7.16, 'O': 5.803,
                      'Si': 4.1491, 'Na': 3.63, 'Cl': 9.577}
scale = 1.0
qpts = np.array(qpts)

def parameter_generator():
    for tt in temp:
        for mat in materials:
            for dw in dw_grid:
                for bos in bose:
                    for ne in negative_e:
                        for cmat in conversion_mat:
                            for lm in lim:
                                yield {'temperature': tt,
                                       'material': mat,
                                       'debye_waller_grid': dw,
                                       'bose': bos,
                                       'negative_e': ne,
                                       'conversion_mat': cmat,
                                       'lim': lm}

def get_expected_output_filename(material_name, pars, opts):
    fname = f"{material_name}_T{pars[0]}"
    if opts.get('debye_waller_grid', None) is not None:
        fname += "_dw" + "".join([str(v) for v in opts['debye_waller_grid']])
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
        if opts['lim']%1 == 0:
            strlim = str(int(opts['lim']))
        else:
            strlim = str(opts['lim'])
        fname += '_lim' + strlim
    fname = fname.replace('.', 'p').replace('-', 'm')
    fname += '.mat'
    return get_abspath(fname, 'expected_output');

def sum_degenerate_modes(w, sf):
    tol = 0.1;
    rows, cols = np.shape(w);
    diff = np.zeros(cols-1);
    summed_sf = np.zeros(np.shape(sf));
    for i in range(rows):
        sum_at = np.where(np.diff(w[i,:]) > tol)[0];
        x = np.zeros(cols);
        x[sum_at + 1] = 1;
        degenerate_modes = np.asarray(np.cumsum(x), dtype=int);
        summed = np.bincount(degenerate_modes, weights=sf[i,:]);
        summed_sf[i, :len(summed)] = summed;
    return summed_sf

def calculate_w_sf(material_pars, material_constructor, par_dict):
    fc = material_constructor(**material_pars)
    par_dict['asr'] = 'reciprocal'
    par_dict['scattering_lengths'] = scattering_lengths
    coherent_sqw = euphonic_sqw_models.CoherentCrystal(fc, **par_dict)
    w, sf = coherent_sqw.horace_disp(qpts[:,0], qpts[:,1], qpts[:,2], scale)
    w = np.array(w).T
    sf = np.array(sf).T
    # Ignore acoustic structure factors by setting to zero - their
    # values can be unstable at small frequencies
    n = int(np.shape(sf)[1] / 2)
    sf[:,:3] = 0.
    if 'negative_e' in par_dict and par_dict['negative_e']:
        sf[:,n:(n+3)] = 0
    # Need to sum over degenerate modes to compare structure factors
    sf_summed = sum_degenerate_modes(w, sf)
    return w, sf_summed

def get_expected_w_sf(fname):
    ref_dat = scipy.io.loadmat(fname)
    expected_w = np.concatenate(ref_dat['expected_w'][0], axis=1)
    expected_sf = np.concatenate(ref_dat['expected_sf'][0], axis=1)
    # Set acoustic structure factors to zero
    expected_sf[:,:3] = 0
    if 'negative_e' in fname:
        n = int(np.shape(expected_sf)[1] / 2)
        expected_sf[:,n:(n+3)] = 0
    expected_sf_summed = sum_degenerate_modes(expected_w, expected_sf)
    return expected_w, expected_sf_summed

@pytest.mark.parametrize("par_dict", parameter_generator())
def test_euphonic_sqw_models(par_dict):
    material_name, material_constructor, material_pars = tuple(
        par_dict.pop('material'))
    expected_w, expected_sf = get_expected_w_sf(
        get_expected_output_filename(
            material_name, [par_dict['temperature'], scale], par_dict))

    for run_par in run_pars:
        par_dict.update(run_par)
        for remove_if_none in list(par_dict.keys()):
            if remove_if_none in par_dict and par_dict[remove_if_none] is None:
                par_dict.pop(remove_if_none)
        w, sf = calculate_w_sf(material_pars, material_constructor, par_dict)
        npt.assert_allclose(w, expected_w, rtol=1e-5, atol=1e-2)
        npt.assert_allclose(sf, expected_sf, rtol=1e-2, atol=1e-2)

# For this test only test a subset of arguments
def scale_test_parameter_generator():
    for tt in temp:
        for mat in materials:
            yield {'temperature': tt,
                   'material': mat,
                   'debye_waller_grid':[6, 6, 6],
                   'negative_e': True,
                   'conversion_mat': (1./2)*np.array([[-1, 1, 1],
                                                      [1, -1, 1],
                                                      [1, 1, -1]])}

# Units of sf have changed, test the calculated and old expected
# values are the same apart from a scale factor
@pytest.mark.parametrize("par_dict", scale_test_parameter_generator())
def test_sf_unit_change(par_dict):
    material_name, material_constructor, material_pars = tuple(
        par_dict.pop('material'))
    fname = get_expected_output_filename(
            material_name, [par_dict['temperature'], scale], par_dict)
    expected_w, expected_sf = get_expected_w_sf(os.path.join(
        os.path.dirname(fname),
        f'old_units_{os.path.basename(fname)}'))

    w, sf = calculate_w_sf(material_pars, material_constructor, par_dict)

    npt.assert_allclose(w, expected_w, rtol=1e-5, atol=1e-2)

    # Change in units angstrom**2 -> mbarn = 1e11. Change from sf per
    # unit cell to sf per atom = 1/n_atoms. Change from relative to
    # absolute units means 1/2 factor is now included
    n_atoms = w.shape[1]/6 # n_atoms is n_modes/6 because of 'negative_e'
    unit_scale = 1e11/(2*n_atoms)
    npt.assert_allclose(sf, expected_sf*unit_scale, rtol=1e-2, atol=1e-2)
