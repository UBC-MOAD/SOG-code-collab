"""Test suite for SOG buildbot diff checker.

Use `python test_check_diff.py` to run the test suite.

:Author: Doug Latornell <dlatorne@eos.ubc.ca>
:Created: 2009-02-02
"""
from __future__ import with_statement
import os
import sys
import StringIO
from tempfile import mkstemp
import unittest
from check_diff import check_diff


class Test_check_diff(unittest.TestCase):
    def assertRaises(self, excClass, callableObj, *args, **kwargs):
        """Monkey-patch unittest.assertRaises() to return exception
        instance.
        """
        try:
            callableObj(*args, **kwargs)
        except excClass, exc_instance:
            return exc_instance
        else:
            if hasattr(excClass,'__name__'): excName = excClass.__name__
            else: excName = str(excClass)
            raise self.failureException, "%s not raised" % excName


    def setUp(self):
        """Pre-test setup.
        """
        # Capture stdout
        sys.stdout = StringIO.StringIO()
        # Create temporary output file stub to use as file1.
        #
        # Note: we have to manage the temporary files ourselves
        # (instead of using tempfile.NamedTemporaryFile) so that the
        # subprocess that runs diff can see them (at least on OS/X).
        handle, self.file1 = mkstemp()
        with open(self.file1, 'w') as fp:
            fp.write("""\
! SOG code standard time series output from physics model
! Time series of iteration count for each time step; mixing layer depth;
! velocity components, temperature, & salinity; at surface, and averaged over
! top 3 m of water column
and surface par
*FromCode: $Source$
*RunDateTime: 2009-03-01 08:29:42
*InitialCTDDateTime: 2003-10-09 11:47:00
*FieldNames: time, iteration count, mixing layer depth, surface u velocity, 3 m avg u velocity, surface v velocity, 3 m avg v velocity, surface temperature, 3 m avg temperature, surface salinity, 3 m avg salinity
surface PAR
""")


    def tearDown(self):
        """Post-test tear-down.
        """
        os.remove(self.file1)
        # Restore stdout to default
        sys.stdout = sys.__stdout__


    def test_too_few_args(self):
        """print usage and exit with status=1 if too few arguments.
        """
        argv = ['check_diff.py']
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), 
                         'usage: %(prog)s file1 file2\n' % {'prog':argv[0]})


    def test_too_many_args(self):
        """print usage and exit with status=1 if too many arguments.
        """
        argv = ['check_diff.py', 'file1', 'file2', 'extra']
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), 
                         'usage: %(prog)s file1 file2\n' % {'prog':argv[0]})


    def test_files_match(self):
        """print nothing and exit with status=0 if files match.
        """
        argv = ['check_diff.py', self.file1, self.file1]
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), '')


    def test_files_different(self):
        """print diff & exit with status=1 if files differ other than date.
        """
        # Create temporary test file with 1 more line than self.file1
        handle, file2 = mkstemp()
        with open(file2, 'w') as fp:
            fp.write("""\
! SOG code standard time series output from physics model
! Time series of iteration count for each time step; mixing layer depth;
! velocity components, temperature, & salinity; at surface, and averaged over
! top 3 m of water column
and surface par
*FromCode: $Source$
*RunDateTime: 2009-03-01 12:29:42
*InitialCTDDateTime: 2003-10-09 11:47:00
*FieldNames: time, iteration count, mixing layer depth, surface u velocity, 3 m avg u velocity, surface v velocity, 3 m avg v velocity, surface temperature, 3 m avg temperature, surface salinity, 3 m avg salinity
surface PAR
spam
""")
        argv = ['check_diff.py', self.file1, file2]
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 1)
        self.assertEqual(sys.stdout.getvalue(), """\
7c7
< *RunDateTime: 2009-03-01 08:29:42
---
> *RunDateTime: 2009-03-01 12:29:42
10a11
> spam

""")
        os.remove(file2)


    def test_files_with_different_RunDateTime(self):
        """print nothing & exit w/ status=0 if files differ only by RunDateTime.
        """
        # Create temporary test file with different RunDateTime line
        # from self.file1
        handle, file2 = mkstemp()
        with open(file2, 'w') as fp:
            fp.write("""\
! SOG code standard time series output from physics model
! Time series of iteration count for each time step; mixing layer depth;
! velocity components, temperature, & salinity; at surface, and averaged over
! top 3 m of water column
and surface par
*FromCode: $Source$
*RunDateTime: 2009-03-01 12:29:42
*InitialCTDDateTime: 2003-10-09 11:47:00
*FieldNames: time, iteration count, mixing layer depth, surface u velocity, 3 m avg u velocity, surface v velocity, 3 m avg v velocity, surface temperature, 3 m avg temperature, surface salinity, 3 m avg salinity
surface PAR
""")
        argv = ['check_diff.py', self.file1, file2]
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), '')
        os.remove(file2)


    def test_files_with_different_Date(self):
        """print nothing & exit w/ status=0 if files differ only by Date.
        """
        # Create temporary file stubs like top of stdout from SOG with
        # different lines
        handle, stdout1 = mkstemp()
        with open(stdout1, 'w') as fp:
            # TODO: capture stdout from SOG
            fp.write("""\
                                              Date = 2009-03-03 02:00:46
                      depth of modelled domain [m] = 40.000000       
                             number of grid points = 80              
                            grid spacing parameter = 0.0000000       
              initialization CTD profile date/time = 2003-10-09 11:47:00 year-day = 282 day-sec = 42420
                              end of run date/time = 2004-11-01 12:00:00 year-day = 306 day-sec = 43200
                                     time step [s] = 900             
                    implicit solver max iterations = 30              
""")
        handle, stdout2 = mkstemp()
        with open(stdout2, 'w') as fp:
            # TODO: change Date line in stdout from SOG
            fp.write("""\
                                              Date = 2009-03-01 17:36:16
                      depth of modelled domain [m] = 40.000000       
                             number of grid points = 80              
                            grid spacing parameter = 0.0000000       
              initialization CTD profile date/time = 2003-10-09 11:47:00 year-day = 282 day-sec = 42420
                              end of run date/time = 2004-11-01 12:00:00 year-day = 306 day-sec = 43200
                                     time step [s] = 900             
                    implicit solver max iterations = 30              
""")
        argv = ['check_diff.py', stdout1, stdout2]
        exception = self.assertRaises(SystemExit, check_diff, argv)
        self.assertEqual(exception.code, 0)
        self.assertEqual(sys.stdout.getvalue(), '')
        os.remove(stdout1)
        os.remove(stdout2)


if __name__ == '__main__':
    unittest.main()
