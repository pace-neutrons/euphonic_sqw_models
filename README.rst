===================
euphonic_sqw_models
===================

A driver for the `Horace <https://github.com/pace-neutrons/Horace.git>`_
inelastic neutron scattering (INS) data analysis package to generate model INS spectra
using the `Euphonic <https://github.com/pace-neutrons/Euphonic.git>`_ package.

This is a pure Python package and is meant to be used in conjunction with
the `horace-euphonic-interface <https://github.com/pace-neutrons/horace-euphonic-interface.git>`_
Matlab package (because Horace is a Matlab program).

`Documentation for the Matlab code <https://horace-euphonic-interface.readthedocs.io/en/latest/>`_.


For developers
==============

Test Data
---------

The test data for both this repository and `horace-euphonic-interface` is stored here,
(since this repo is included as a submodule in `horace-euphonic-interface`)
but in Matlab `.mat` file format.
This means that Matlab should be used to generate it if it becomes outdated.
Use `runtests('test/EuphonicGenerateTestData.m')` in the `horace-euphonic-interface` folder.
(You may also have to set `generate_test_data` to `true` in that file.
Then copy the `*.mat` files from `test/expected_output` to the same folder in this repository
and commit.


Euphonic version
----------------

The minimum Euphonic version is set in the `min_requirements.txt` file.
This information will be inherited by `horace-euphonic-interface`.
