import unittest
from fs_helpers import *

TEST_FILE_PATH = "./test_data/test_file.txt"
NOT_EXISTS_FILE_PATH = "not_exists.txt"


class ReadFileTest(unittest.TestCase):

    def test_read_file_read_content(self):
        """ read_file function test """
        content = read_file(TEST_FILE_PATH)
        self.assertEqual(content, "test1\ntest2", "file content read correctly")

    def test_read_file_file_not_exists(self):
        try:
            read_file(NOT_EXISTS_FILE_PATH)
        except FileNotFoundError:
            pass
        else:
            fail("expected a FileNotFoundError")

    def test_read_lines(self):
        """ read_lines function test """
        lines = read_lines(TEST_FILE_PATH)
        self.assertEqual(lines, ["test1", "test2"], "lines read read correctly")

    def test_read_lines_file_not_found(self):
        try:
            read_lines(NOT_EXISTS_FILE_PATH)
        except FileNotFoundError:
            pass
        else:
            fail("expected a FileNotFoundError")


if __name__ == "__main__":
    unittest.main()
