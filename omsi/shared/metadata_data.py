"""
Define a dict for describing metadata
"""

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
        dict.__setitem__(self, key, value)

    def get_metadata_values(self):
        """
        Get a list of all metadata values.
        """
        return [metadata_value['value'] for metadata_value in self.values()]

    def get_metadata_descriptions(self):
        """
        Get a list of all metadata descriptions
        """
        return [metadata_value['description'] for metadata_value in self.values()]

    def get_metadata_units(self):
        """
        Get a list of all metadata units
        """
        return [metadata_value['unit'] for metadata_value in self.values()]


class metadata_value(dict):
    """
    A single metadata value
    """
    def __init__(self, value, description, unit=None):
        super(metadata_value, self).__init__()
        dict.__setitem__(self, 'value', value)
        dict.__setitem__(self, 'description', description)
        dict.__setitem__(self, 'unit', unit)

    def __setitem__(self, key, value):
        if key in self.keys():
            dict.__setitem__(self, key, value)
        else:
            raise KeyError('Illegal key ' + str(key) + ' given. Valid keys are ' + str(self.keys()))

