import numpy


class Relation(object):
    """A Relation object has a pair of NumPy arrays containing the
    independent and dependent data values of a data set.  It also has
    attributes that contain the filespec from which the data is read,
    and the units of the data arrays.

    This is a base class for implementing specific types of relations
    such as timeseries, or profiles.  This class provides a
    read_data() method, but expects classes that inherit from it to
    implement a _read_header() method that returns lists of field
    names and field units read from the SOG data file header.
    """
    def __init__(self, datafile):
        """Create a Relation instance with its datafile attribute initialied.
        """
        self.datafile = datafile


    def _read_header(self):
        """This method is expected to be implemented by classes based
        on Relation.

        It must return lists of field names and field units from the
        SOG data file header.
        """
        raise NotImplementedError

    
    def read_data(self, indep_field, dep_field):
        """Read the data for the specified independent and dependent
        fields from the data file.

        Sets the indep_data and dep_data attributes to NumPy arrays,
        and the indep_units and dep_units attributes to units strings
        for the data fields.
        """
        try:
            f = open(self.datafile, 'r')
        except IOError, msg:
            return_code = 2
            print msg
            raise SystemExit(return_code)
        # Read the file header
        (field_names, field_units) = self._read_header(f)
        # Translate the field names into column numbers
        for field in (indep_field, dep_field):
            if field not in field_names:
                return_code = 2
                msg = ("Invalid field name: %(field)s, "
                       "choices are:\n %(field_names)s" % vars())
                print msg
                raise SystemExit(return_code)
        indep_col = field_names.index(indep_field)
        dep_col = field_names.index(dep_field)
        # Data units
        self.indep_units = field_units[indep_col]
        self.dep_units = field_units[dep_col]
        # Read the data
        self.indep_data, self.dep_data = [], []
        try:
            for line in f:
                self.indep_data.append(float(line.split()[indep_col]))
                self.dep_data.append(float(line.split()[dep_col]))
        finally:
            f.close()
        # Transform data lists into NumPy arrays
        self.indep_data = numpy.array(self.indep_data)
        self.dep_data = numpy.array(self.dep_data)
