`Unreleased <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v1.0.0...HEAD>`_
----------

 - New Features:

   - ``CoherentCrystal`` can now also be created with a
     ``euphonic.brille.BrillerInterpolator`` object

`v1.0.0 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.5.2...v1.0.0>`_
------

- Changes:

  - ``n_threads`` is now explicitly named in the ``CoherentCrystal`` constructor arguments
    rather than being part of ``**kwargs``
  - ``psutil`` has been added as a dependency

- Improvements:

  - If ``chunk`` hasn't been provided to ``CoherentCrystal``, a recommended chunk
    size will now be calculated and set based on available memory.

`v0.5.2 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.5.1...v0.5.2>`_
------

- Improvements:

  - Improve reliability by testing on oldest and newest supported Euphonic versions

`v0.5.1 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.5.0...v0.5.1>`_
----------

- Bug fixes:

  - Use of temperature=0 will now calculate the 0K Debye-Waller and Bose
    population factors - previously these temperature dependent effects
    were not calculated at 0K

`v0.5.0 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.4.0...v0.5.0>`_
------

- Improvements:

  - There is a new ``frequency_scale=1.0`` argument to ``horace_disp`` which
    allows the output frequencies to be scaled

- Breaking changes:

  - The ``pars=[]`` argument to ``horace_disp`` has been changed to
    ``intensity_scale=1.0``

`v0.4.0 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.3.0...0.4.0>`_
------

- Dependency changes:

  - Euphonic version dependency increased from >= 0.5.0 to >= 0.6.0

- Breaking changes:

  - The default units of ``StructureFactor.structure_factors`` in Euphonic have been
    changed from ``angstrom**2`` per unit cell to ``mbarn`` per sample atom, and are
    now in absolute units including a previously omitted 1/2 factor. So the structure
    factors produced by CoherentCrystal.horace_disp have increased by a factor of
    ``1e11/(2*N_atoms)``

- Other changes:

  - The ``eta_scale`` keyword argument to ``CoherentCrystal`` has been deprecated,
    ``dipole_parameter`` should be used instead
  - A ValueError will now be raised if an unrecognised keyword argument is passed
    to ``CoherentCrystal``


`v0.3.0 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.2...v0.3.0>`_
------

- Dependency changes:

  - Euphonic version dependency increased from >=0.4.0 to >=0.5.0

- Breaking changes:

  - ``fall_back_on_python`` argument to ``horace_disp`` has been removed as this has
    been removed in Euphonic

`v0.2 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/v0.1.0...v0.2>`_
------

- Breaking changes:

  - Major update to how Horace-Euphonic-Interface works, most code has been rewritten in
    Python to allow it to be called directly from the Python version of Horace. It has also
    been split into two separate repositories so that the Python version of Horace only
    needs to include what it needs, and to allow easier updating and management of
    either Python or Matlab parts of the code. The two repositories are:

     - This repository (``euphonic_sqw_models``), which contains the Python part of the code
     - `Horace-Euphonic-Interface <https://github.com/pace-neutrons/horace-euphonic-interface>`_,
       which has retained its name, but now only includes minimal Matlab wrappers around
       the Python code in this repository.

   - There has also been a major refactor, the main changes are:

     - ``euphonic_sf`` has been removed
     - ``euphonic_on`` has been removed
     - Force constants are now a separate object (``ForceConstants``) rather than
       passing these arguments to ``euphonic_sf``
     - The model parameters are set in a ``CoherentCrystal`` model object, rather than
       passing these parameters to ``euphonic_sf``
     - The function handle to be passed to ``disp2sqw_eval`` is ``CoherentCrystal.horace_disp``
       rather than ``euphonic_sf``
     - The ``dw_grid`` argument has been renamed to ``debye_waller_grid``

- Other changes:

  - A ``verbose`` argument has been added to ``horace_disp`` which can be set to ``Flase``
    to prevent printing of q-point progress


`v0.1.0 <https://github.com/pace-neutrons/euphonic_sqw_models/compare/81607231b...v0.1.0>`_
------

- First release
