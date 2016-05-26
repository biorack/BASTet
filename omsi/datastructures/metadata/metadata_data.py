"""
Define infrastructure for describing metadata (in memory)
"""
import json

class metadata_dict(dict):
    """
    Dictionary's for storing metadata information. The values in the dict must be
    of type omsi_metadata_value and the keys must be strings.
    """
    def __init__(self):
        super(metadata_dict, self).__init__()

    def __setitem__(self, key, value):
        if not isinstance(key, basestring):
            raise KeyError('Metadata keys must be string names')
        if not isinstance(value, metadata_value):
            raise ValueError('Metadata values must be defined via metadata_value dicts')
        if value['name'] is None:
            value['name'] = key
        dict.__setitem__(self, key, value)

    def get_metadata_values(self):
        """
        Get a list of all metadata values.
        """
        return [metadata_val['value'] for metadata_val in self.values()]

    def get_metadata_descriptions(self):
        """
        Get a list of all metadata descriptions
        """
        return [metadata_val['description'] for metadata_val in self.values()]

    def get_metadata_units(self):
        """
        Get a list of all metadata units
        """
        return [metadata_val['unit'] for metadata_val in self.values()]


class metadata_value(dict):
    """
    A single metadata value

    Allowed keys:
    -------------

    * `name` The name of the metadata value
    * `value` The actual value associated with the metadata object
    * `description` The text description of the metadata object
    * `unit` The unit string
    * `ontology` Optional ontology

    Data values:
    ------------
    - If data values are dicts, then the dict should  be JSON serializable in order for
      the value to be stored in omsi HDF5 files
    - If data values are lists, then the lists must be able to be converted to numpy
      in a way that h5py understands


    """
    def __init__(self,
                 name,
                 value,
                 description,
                 unit=None,
                 ontology=None):
        """

        :param name: The name of the metadata value. Name may be None if the metadata_value
            is added to a metadata_dict as it will be set (if missing) when adding it to the
            metadata_dict
        :param value: The actual value associated with the metadata object
        :param description: The text description of the metadata object
        :param unit: The unit string
        :param ontology: Optional ontology

        """
        super(metadata_value, self).__init__()
        # Initalize all the keys
        dict.__setitem__(self, 'value', None)
        dict.__setitem__(self, 'description', None)
        dict.__setitem__(self, 'unit', None)
        dict.__setitem__(self, 'ontology', None)
        dict.__setitem__(self, 'name', None)

        # Set all the values for the keys
        self['value'] = value
        self['description'] = description
        self['unit'] = unit
        self['ontology'] = ontology
        self['name'] = name

    def __setitem__(self, key, value):
        if key in self.keys():
            if key == 'name' and not (isinstance(value, basestring) or value is None):
                raise ValueError('The name value must be a string')
            elif key == 'description' and not (isinstance(value, basestring) or value is None):
                raise ValueError('The description value must be a string')
            elif key == 'unit' and not (isinstance(value, basestring) or value is None):
                raise ValueError('The unit must be a string')

            dict.__setitem__(self, key, value)
        else:
            raise KeyError('Illegal key ' + str(key) + ' given. Valid keys are ' + str(self.keys()))

    def __getattr__(self, item):
        """
        Expose the keys also as attributes for usability.
        :param item: Name of the attribute
        :return: Value of the attribute
        :raises: AttributeError if the attribute is not found
        """
        if item in self.keys():
            return self[item]
        else:
            raise AttributeError("'metadata_value' has not attribute " + str(item))

