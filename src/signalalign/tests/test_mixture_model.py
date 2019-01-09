#!/usr/bin/env python
"""Testing mixture_model_test.py"""
########################################################################
# File: mixture_model_test.py
#
# Author: Andrew Bailey
# History: 1/7/19 Created
########################################################################


import os
import unittest


class SeqToolsTests(unittest.TestCase):
    """Test the functions in seq_tools.py"""

    @classmethod
    def setUpClass(cls):
        super(SeqToolsTests, cls).setUpClass()
        cls.HOME = '/'.join(os.path.abspath(__file__).split("/")[:-1])

    def test_something(self):
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
