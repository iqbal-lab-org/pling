import unittest

# Instance of the TestCase to use the assertEqual method
tc = unittest.TestCase()


def read_file(path):
    """Helper method to read a file."""
    with open(path, 'rb') as f:
        return f.read()


def assert_files_are_identical(path_to_first_file, path_to_second_file):
    first_file_content = read_file(path_to_first_file)
    second_file_content = read_file(path_to_second_file)
    tc.assertEqual(first_file_content, second_file_content)