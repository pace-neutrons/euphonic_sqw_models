import warnings
from typing import Optional, Dict, Union, Tuple

import psutil
import numpy as np

from euphonic import (ForceConstants, QpointPhononModes, DebyeWaller,
                      Quantity, ureg)
from euphonic.util import mp_grid

# Allow Euphonic's deprecation warnings to be raised if deprecated
# arguments are passed through this module
warnings.filterwarnings('default', category=DeprecationWarning,
                        module=__name__)

class CoherentCrystal(object):
    """
    A class to calculate neutron scattering dispersions surfaces for
    use in Horace. It must be constructed from a Euphonic
    ForceConstants or euphonic.brille.BrilleInterpolator object.

    Attributes
    ----------
    force_constants : Union[ForceConstants, BrilleInterpolator]
        An object that has a calculate_qpoint_phonon_modes function
        that produces a euphonic.QpointPhononModes object. Typically
        ForceConstants or BrilleInterpolator.
    debye_waller_grid : (3,) float ndarray or None (default: None)
        The Monkhorst-Pack q-point grid size to use for the
        Debye-Waller calculation
    debye_waller : DebyeWaller object
        This is computed and set internally by the class when the
        computation requests it (debye_waller_grid is set and
        temperature is not `None`)
    temperature : float, float Quantity or None (default: None)
        Temperature at which to calculate phonons (used for Bose and
        Debye-Waller factor calculations). If None no temperature
        dependent effects will be applied
    bose : bool (default: True)
        Whether to include the Bose temperature factor in the
        calculated structure factor
    negative_e : bool (default: False)
        Whether to calculate the phonon annihilation / neutron energy
        gain / negative energy spectrum
    conversion_mat : (3, 3) float ndarray or None (default: None)
        This matrix is applied to the input (hkl) q-point values before
        the phonon calculation
    chunk : int (default: None)
        Will chunk the phonon calculation into blocks of
        <chunk> q-points. If not provided, an optimum chunk
        size will be calculated based on available memory. If
        a chunk size of zero is provided, the calculation will
        not be chunked.
    lim : float (default: infinity)
        A cut-off for the calculated structure factor where values
        above `lim` are set equal to it
    scattering_lengths : dict or string (default: 'Sears1992')
        The scattering lengths to use for the phonon calculations
        either as a dictionary with elements as keys and floats
        (in fm) or Quantities as values, or a string denoting a
        database (internal or from a file)
    verbose : bool (default: True)
        Whether to print information on calculation progress
    calc_modes_kwargs
        Any other keyword arguments (e.g asr, n_threads) will be
        passed to euphonic.ForceConstants.calculate_qpoint_phonon_modes.
        See the available arguments with
        help(euphonic.ForceConstants.calculate_qpoint_phonon_modes)
    """
    def __init__(
            self,
            force_constants: Union[ForceConstants, 'BrilleInterpolator'],
            debye_waller_grid: Optional[Tuple[int, int, int]] = None,
            debye_waller: Optional[DebyeWaller] = None,
            temperature: Union[float, Quantity, None] = None,
            bose: bool = True,
            negative_e: bool = False,
            conversion_mat: Optional[np.ndarray] = None,
            chunk: int = None,
            lim: float = np.inf,
            scattering_lengths: Union[
                str, Dict[str, Union[Quantity, float]]] = 'Sears1992',
            verbose: bool = True,
            **calc_modes_kwargs) -> None:
        """
        Parameters
        ----------
        force_constants
            An object that has a calculate_qpoint_phonon_modes function
            that produces a euphonic.QpointPhononModes object.
            Typically ForceConstants or BrilleInterpolator.
        debye_waller_grid
            The Monkhorst-Pack q-point grid size to use for the
            Debye-Waller calculation. If None the Debye-Waller
            calculation is not performed.
        debye_waller
            This is computed and set internally by the class when the
            computation requests it (debye_waller_grid is set and
            temperature is not None)
        temperature
            Temperature at which to calculate phonons (used for Bose and
            Debye-Waller factor calculations). If None no
            temperature-dependent effects will be applied
        bose
            Whether to include the Bose temperature factor in the
            calculated structure factor
        negative_e
            Whether to calculate the phonon annihilation / neutron energy
            gain / negative energy spectrum
        conversion_mat
            Shape (3, 3) float ndarray. This matrix is applied to the
            input (hkl) q-point values before the phonon calculation
        chunk : int (default: None)
            Will chunk the phonon calculation into blocks of
            <chunk> q-points. If not provided, an optimum chunk
            size will be calculated based on available memory. If
            a chunk size of zero is provided, the calculation will
            not be chunked.
        lim
            A cut-off for the calculated structure factor where values
            above `lim` are set equal to it
        scattering_lengths
            The scattering lengths to use for the phonon
            calculations either as a dictionary with elements as keys and
            floats (in fm) or Quantities as values, or a
            string denoting a database (internal or from a file)
        verbose
            Whether to print information on calculation progress
        **calc_modes_kwargs
            Any other keyword arguments (e.g asr, use_c) will be
            passed to the calculate_qpoint_phonon_modes method of
            the input ForceConstants or BrilleInterpolator object.
        """
        self.force_constants = force_constants
        self.debye_waller_grid = debye_waller_grid
        self.debye_waller = debye_waller
        self.temperature = temperature
        self.bose = bose
        self.negative_e = negative_e
        self.conversion_mat = conversion_mat
        if chunk is None:
            chunk = self.get_optimum_chunk_size(
                force_constants.crystal.n_atoms)
        self.chunk = chunk
        self.lim = lim
        self.scattering_lengths = scattering_lengths
        self.verbose = verbose
        # Explicitly convert to integers so users calling from MATLAB
        # don't have to
        int_kwargs = ['threads', 'n_threads']
        for ikw in int_kwargs:
            if ikw in calc_modes_kwargs:
                calc_modes_kwargs[ikw] = int(calc_modes_kwargs[ikw])
        self.calc_modes_kwargs = calc_modes_kwargs

    def _calculate_sf(self, qpts: np.ndarray
                      ) -> Tuple[np.ndarray, np.ndarray]:
        if self.temperature is not None:
            if (self.debye_waller is None
                and self.debye_waller_grid is not None):
                self._calculate_debye_waller()
        phonons = self._calculate_phonon_modes(qpts)
        sf_obj = phonons.calculate_structure_factor(
            scattering_lengths=self.scattering_lengths,
            dw=self.debye_waller)
        w = sf_obj.frequencies.magnitude
        sf = sf_obj.structure_factors.magnitude
        if self.temperature is not None and self.bose:
            bose = sf_obj._bose_factor(self.temperature)
            sf = (1 + bose) * sf
        if self.negative_e:
            w = np.hstack((w, -w))
            neg_sf = sf_obj.structure_factors.magnitude
            if self.temperature is not None and self.bose:
                neg_sf = bose * neg_sf
            sf = np.hstack((sf, neg_sf))
        return w, sf

    def horace_disp(self, qh: np.ndarray, qk: np.ndarray, ql: np.ndarray,
                    intensity_scale: float = 1.0,
                    frequency_scale: float = 1.0
                    ) -> Tuple[Tuple[np.ndarray, ...],
                               Tuple[np.ndarray, ...]]:
        """
        Calculates the phonon dispersion surface for input qh, qk, and
        ql vectors for use with Horace

        Parameters
        ----------
        qh
            The q-points in H to calculate
        qk
            The q-points in K to calculate
        ql
            The q-points in L to calculate
        intensity_scale
            The factor to multiply the intensity by
        frequency_scale
                The factor to multiply the phonon frequencies
                by, as DFT can often over/underestimate frequencies

        Returns
        -------
        w
            Length (n_qpts,) tuple of length (n_modes,) float ndarrays.
            The phonon dispersion energies
        sf
            Length (n_qpts,) tuple of length (n_modes,) float ndarrays.
            The dynamical structure corresponding to phonon energies in
            w
        """
        if self.chunk > 0:
            lqh = len(qh)
            for i in range(int(np.ceil(lqh / self.chunk))):
                qi = i * self.chunk
                qf = min((i+1) * self.chunk, lqh)
                if self.verbose:
                    print(f'Using Euphonic to interpolate for q-points '
                          f'{qi}:{qf} out of {lqh}')
                qpts = np.stack((np.squeeze(qh[qi:qf]),
                                 np.squeeze(qk[qi:qf]),
                                 np.squeeze(ql[qi:qf])),
                                 axis=1)
                if self.conversion_mat is not None:
                    qpts = np.matmul(qpts, self.conversion_mat)
                sqw = self._calculate_sf(qpts)
                if i == 0:
                    w = sqw[0]
                    sf = sqw[1]
                else:
                    w = np.vstack((w, sqw[0]))
                    sf = np.vstack((sf, sqw[1]))
        else:
            qpts = np.stack((np.squeeze(qh),
                             np.squeeze(qk),
                             np.squeeze(ql)),
                             axis=1)
            if self.conversion_mat is not None:
                qpts = np.matmul(qpts, self.conversion_mat)
            w, sf = self._calculate_sf(qpts)
        if frequency_scale != 1.:
            w *= frequency_scale
        if intensity_scale != 1.:
            sf *= intensity_scale
        sf = np.minimum(sf, self.lim)
        # Splits into different dispersion surfaces
        # (python tuple == matlab cell)
        # But the data must be contiguous in memory so we need to do a
        # real tranpose (.T just changes strides)
        # So we need to convert to "fortran" format (which physically
        # transposes data) before doing ".T"
        w = np.asfortranarray(w).T
        sf = np.asfortranarray(sf).T
        return tuple(w), tuple(sf)

    def _calculate_phonon_modes(self, qpts: np.ndarray
                                ) -> QpointPhononModes:
        return self.force_constants.calculate_qpoint_phonon_modes(
            qpts, **self.calc_modes_kwargs)

    def _calculate_debye_waller(self) -> None:
        if self.debye_waller_grid is None:
            raise RuntimeError(
                'Q-points grid for Debye Waller calculation not set')
        dw_qpts = mp_grid(self.debye_waller_grid)
        dw_phonons = self._calculate_phonon_modes(dw_qpts)
        self.debye_waller = dw_phonons.calculate_debye_waller(
            self.temperature)

    @staticmethod
    def get_optimum_chunk_size(n_atoms: int):
        mem = psutil.virtual_memory().available
        evec_bytes_per_qpt = 16*(3*n_atoms)**2
        # Divide chunk by 10 for a conservative estimate
        chunk = int(mem/(10*evec_bytes_per_qpt))
        return chunk

    @property
    def force_constants(self) -> Union[ForceConstants, 'BrilleInterpolator']:
        return self._force_constants

    @force_constants.setter
    def force_constants(self, val: Union[ForceConstants, 'BrilleInterpolator']
                        ) -> None:
        if hasattr(val, 'calculate_qpoint_phonon_modes'):
            self._force_constants = val
        else:
            raise ValueError('Invalid force constant model')

    @property
    def debye_waller(self) -> Union[DebyeWaller, None]:
        return self._debye_waller

    @debye_waller.setter
    def debye_waller(self, val: Union[DebyeWaller, None, str]) -> None:
        if val is None or val == 'None':
            self._debye_waller = None
        else:
            if hasattr(val, 'debye_waller'):
                self._debye_waller = val
            else:
                raise ValueError('Invalid Debye-Waller object')

    @property
    def debye_waller_grid(self) -> Union[Tuple[int, int, int], None]:
        return self._debye_waller_grid

    @debye_waller_grid.setter
    def debye_waller_grid(self, val: Union[Tuple[int, int, int], None, str]
                          ) -> None:
        if val is None or val == 'None':
            self._debye_waller_grid = None
        else:
            val = tuple(np.squeeze(np.array(val)))
            if np.shape(val) == (3,):
                self._debye_waller_grid = tuple([int(v) for v in val])
                # Reset the Debye Waller factor if it was previously set.
                self.debye_waller = None
            else:
                raise ValueError('Invalid Debye-Waller grid')

    @property
    def conversion_mat(self) -> Union[np.ndarray, None]:
        return self._conversion_mat

    @conversion_mat.setter
    def conversion_mat(self, val: Union[np.ndarray, None, str]) -> None:
        if val is None or (isinstance(val, str) and val.startswith('None')):
            self._conversion_mat = None
        else:
            val = np.array(val)
            if np.shape(val) == (3, 3):
                self._conversion_mat = val
            else:
                raise ValueError('Invalid conversion matrix')

    @property
    def temperature(self) -> Union[Quantity, None]:
        return self._temperature

    @temperature.setter
    def temperature(self, val: Union[Quantity, float, None]) -> None:
        if val is None or val == 'None':
            self._temperature = None
        else:
            if isinstance(val, (float, int)):
                val = val*ureg('K')
            if hasattr(val, 'units') and val >= 0*ureg('K'):
                self._temperature = val
            else:
                raise ValueError('Invalid temperature')

    @property
    def scattering_lengths(self) -> Union[str, Dict[str, Quantity]]:
        return self._scattering_lengths

    @scattering_lengths.setter
    def scattering_lengths(
            self, val: Union[str, Dict[str, Union[float, Quantity]]]) -> None:
        if isinstance(val, str):
            self._scattering_lengths = val
        elif isinstance(val, dict):
            self._scattering_lengths = \
                {ky: (v if hasattr(v, 'units') else v * ureg('fm'))
                     for ky, v in val.items()}
        else:
            raise ValueError('Invalid scattering lengths')

    @property
    def chunk(self) -> int:
        return self._chunk

    @chunk.setter
    def chunk(self, val: int) -> None:
        self._chunk = int(val)
