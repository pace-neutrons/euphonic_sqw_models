import pytest
import coverage
import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--coverage", dest="coverage",
        help="Report coverage in a file with this name")
    parser.add_argument(
        "-r", "--report", action="store_true",
        help="If present report test results to junit_report*.xml files")
    parser.add_argument(
        "-t", "--test-file", dest="test_file", action="store",
        help="The test file to run", default="test")
    parser.add_argument(
        "-m",
        help="Run tests with this marker (is passed directly to pytest)")
    args_parsed = parser.parse_args()

    pytest_options = [args_parsed.test_file]
    if args_parsed.m:
        pytest_options += ['-m', args_parsed.m]
    if args_parsed.report:
        pytest_options.append('--junitxml=report.xml')
    if args_parsed.coverage:
        cov = coverage.Coverage(source=['euphonic_sqw_models'], omit=['*/_version.py'])
        cov.start()
    pytest_options.insert(0, '-v')
    test_output = pytest.main(pytest_options)
    
    if args_parsed.coverage:
        cov.stop()
        cov.xml_report(outfile=args_parsed.coverage)

    sys.exit(test_output)
