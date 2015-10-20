"""
Define ontologies for metadata

:var METADATA_ONTOLOGIES: Description of simple ontologies for metadata. This is a dict where the
     values are the descriptions of the ontologies and the key is the commonly used name of the
     variable associated with the ontology, however, in practice an ontology may be associated
     with many different metadata values (the purpose of th ontology is to standardize the values
     not the names of metadata variables). Available ontologies include:

    * `polarity` : Description of the polarity of the instrument
    * `msn_value_of_n` : Numeric level of mass spectrometry used (e.g., 1 for MS1 etc).

An ontology can be many things. In general an ontology is a mechanism to formally name and define the
types, properties, and interrelationships of entities. We here refer broadly to the concept of ontologies
as a means to standardize the names of things. To define an ontology we use simple dicts with the
following key/value entries:

    * `name`: The name of the ontology
    * `value` : Python dict describing the actual ontology. Often this is simply a dict where the \
                the keys are the allowed values and the value is the textual description of the \
                meaning of that value. In many cases this may be a more complete description of \
                an ontology.
    * `unit` : The standard unit associated with the values (or None if not unit is available)
    * `description` : Human-readable textual description of the ontology
    * `version` : The version of the ontology used
    * 'uri` : The Universal Remote Identifier (often a URL) associated with the ontology (or None)

Ontologies are stored in the OpenMSI data format as JSON attributes associated with the metadata, i.e.,
one main restriction is that ontologies should be JSON serializable (which in most cases should be a
problem).

"""

METADATA_ONTOLOGIES = {}

METADATA_ONTOLOGIES['polarity'] = {'name': 'bo_ms_polarity',
                                   'value': {'pos': 'positive polarity',
                                             'neg': 'negative_polarity'},
                                   'unit': None,
                                   'description': 'Polarity',
                                   'version': '0.1.0',
                                   'URI': None}

METADATA_ONTOLOGIES['msn_value_of_n'] = {'name': 'bo_msn_value_of_n',
                                         'value': {1: 'MS1', 2: 'MS2', 3: 'MS3', 4: 'MS4',
                                                   5: 'MS5', 6: 'MS6', 7: 'MS7', 8: 'MS8',
                                                   9: 'MS9', 10: 'MS10'},
                                         'unit': 'level',
                                         'description': 'The level n of MS^n mass spectrometry',
                                         'version': '0.1.0',
                                         'URI': None}


class metadata_ontologies(dict):
    """
    Helper class for interacting with ontologies
    """
    def __init__(self):
        super(metadata_ontologies, self).__init__()
        self.update(METADATA_ONTOLOGIES)

