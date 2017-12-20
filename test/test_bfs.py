import os
import unittest

import pandas as pd
from pandas.util import testing as pdt

from physep import hysep

pd.set_option('display.precision', 7)
pd.set_option('display.expand_frame_repr', False)
TEST_DATA_FILENAME = os.path.join(os.path.dirname(__file__), 'hysep_data.csv')
KEYWORDS = {'size': 93}


class TestHySepMethods(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestHySepMethods, self).__init__(*args, **kwargs)
        self.test_data = pd.Series.from_csv(TEST_DATA_FILENAME,
                                            parse_dates=True,
                                            index_col=0,
                                            header=0)
        self.filters = [key for key in hysep.__dict__ if key.endswith('filter')]

    def test_input(self):
        self.assertIsInstance(self.test_data, pd.Series,
                              'test data is not type pd.Series')
        self.assertIsInstance(self.test_data.index, pd.DatetimeIndex,
                              'test data index is not type pd.DateTimeIndex')

    def test_minimum(self):
        self.generic_bfs_test(hysep.minimum_filter)

    def test_fixed_interval(self):
        self.generic_bfs_test(hysep.fixed_interval_filter)

    def test_sliding_interval(self):
        self.generic_bfs_test(hysep.sliding_interval_filter)

    def test_local_minimum(self):
        print 'morgane'
        self.generic_bfs_test(hysep.local_minimum_filter)

    def generic_bfs_test(self, bsf_filter):
        baseflow, quickflow = bsf_filter(self.test_data, size=KEYWORDS['size'])

        # check lengths
        len_arr = [len(x) for x in [self.test_data, baseflow, quickflow]]
        self.assertTrue(len_arr.count(len_arr[0]) == len(len_arr))

        # check indexes
        pdt.assert_index_equal(self.test_data.index, baseflow.index)
        pdt.assert_index_equal(baseflow.index, quickflow.index)

        # check series val
        self.assertTrue((baseflow >= 0).any())
        self.assertTrue((quickflow >= 0).any())

        # check series sum
        sumflow = baseflow + quickflow
        df_na = pd.concat([self.test_data, sumflow], axis=1).dropna()

        pdt.assert_series_equal(sumflow, self.test_data,
                                check_dtype=True,
                                check_index_type=True,
                                check_less_precise=15,
                                check_names=False)


if __name__ == '__main__':

    unittest.main()
