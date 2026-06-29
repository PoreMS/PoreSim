import os
import pytest


@pytest.fixture(autouse=True, scope="session")
def change_to_tests_dir():
    """Run all tests from within the tests/ directory so relative paths work."""
    original = os.getcwd()
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(tests_dir)
    yield
    os.chdir(original)
