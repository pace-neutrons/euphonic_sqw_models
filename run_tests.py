import pytest
import coverage
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--coverage", action="store_true",
        help="If present report coverage in a coverage.xml file in reports")
    parser.add_argument(
        "-r", "--report", action="store_true",
        help="If present report test results to junit_report*.xml files")
    parser.add_argument(
        "-t", "--test-file", dest="test_file", action="store",
        help="The test file to run", default="test")
    args_parsed = parser.parse_args()

    pytest_options = [args_parsed.test_file]
    if args_parsed.report:
        pytest_options.append('--junitxml=report.xml')
    
    if args_parsed.coverage:
        cov = coverage.Coverage(source=['euphonic_horace'], omit=['*/_version.py'])
        cov.start()

    test_output = pytest.main(pytest_options)
    
    if args_parsed.coverage:
        cov.stop()
        cov.xml_report(outfile='coverage.xml')
