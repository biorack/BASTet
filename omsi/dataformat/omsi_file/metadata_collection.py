"""
Module for management of general metadata storage entities. These are often specialized
---e.g., omsi_file_instrument, omsi_file_sample---to store specific metadata and add
more functionality.
"""
import time
import json

import h5py
from numpy import nan as float_nan
import numpy as np

from omsi.dataformat.omsi_file.format import omsi_format_common, omsi_format_metadata_collection
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager
from omsi.datastructures.metadata.metadata_data import metadata_dict, metadata_value
from omsi.shared.log import log_helper


class omsi_metadata_collection_manager(omsi_file_object_manager):
    """
    This is a file format manager helper class
    used to define common functionality needed for management of metadata-related data.
    Usually, a class that defines a format that contains an omsi_file_metadata object
    will inherit from this class (in addition to omsi_file_common) to acquire the common
    features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar metadata_parent: The parent h5py.Group object containing the method object to be managed

    """
    def __init__(self, metadata_parent=None):
        super(omsi_metadata_collection_manager, self).__init__(metadata_parent)
        self.metadata_parent = metadata_parent

    def create_metadata_collection(self, group_name=None, metadata=None, flush_io=True):
        """
        Add a new group for managing metadata

        :param group_name: Optional name of the new metadata group. If None is given
            then the omsi_format_metadata_collection.metadata_collection_groupname_default will be used
        :type group_name: str, None
        :param metadata: Additional metadata to be added to the collection after creation
        :type metadata: None, omsi.shared.metadata_data.metadata_value, omsi_shared.metadata_data.metadata_dict
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 buffers are
                         flushed so that all data has been written to file

        :returns: `omsi_file_metadata_collection` object
        """
        return omsi_file_metadata_collection.___create___(self.metadata_parent,
                                                          group_name=group_name,
                                                          metadata=metadata,
                                                          flush_io=flush_io)

    def get_metadata_collections(self, omsi_object=None, name=None):
        """
        Get all metadata_collections defined for given OpenMSI file API object or h5py.Group.

        :param omsi_object: The omsi file API object or h5py.Group object that we should check.
            If set to None (default) then the self.metadata_parent will be use
        :param name: If name is specified, then only retrieve collections with the given name

        :returns: List of omsi_file_metadata_collection objects for the requested group. The function
            returns None in case that the h5py.Group for the omsi_object could not be determined.

        """
        metadata_collections = []
        use_omsi_object = self.metadata_parent if omsi_object is None else omsi_object
        h5py_group = None
        if isinstance(use_omsi_object, omsi_file_common):
            h5py_group = use_omsi_object.managed_group
        elif isinstance(use_omsi_object, h5py.Group):
            h5py_group = use_omsi_object
        elif isinstance(use_omsi_object, h5py.Dataset):
            h5py_group = use_omsi_object.parent
        elif isinstance(use_omsi_object, h5py.File):
            h5py_group = use_omsi_object['/']

        if h5py_group is not None:
            for h5py_name, h5py_value in h5py_group.iteritems():
                curr_omsi_object = omsi_file_common.get_omsi_object(h5py_value)
                if isinstance(curr_omsi_object, omsi_file_metadata_collection):
                    if name is None or h5py_name == name:
                        metadata_collections.append(curr_omsi_object)
        else:
            return None

        return metadata_collections

    def has_metadata_collections(self, omsi_object=None):
        """
        Check whether the given omsi API object (or h5py.Group) contains any metadata collections

        :param omsi_object: The omsi file API object or h5py.Group object that we should check.
            If set to None (default) then the self.metadata_parent will be used

        :returns: Boolean indicating whether metadata collections were found
        """
        collections = self.get_metadata_collections(omsi_object=omsi_object)
        found_metadata_collection = collections is not None
        if found_metadata_collection:
            found_metadata_collection = len(collections) > 0
        return found_metadata_collection

    def get_default_metadata_collection(self, omsi_object=None):
        """
        Get the default metadata collection object if it exists

        :param omsi_object: The omsi file API object or h5py.Group object that we should check.
            If set to None (default) then the self.metadata_parent will be used
        :return: None, omsi_file_metadata_collection
        """
        default_collections = self.get_metadata_collections(
            omsi_object=omsi_object,
            name=omsi_format_metadata_collection.metadata_collection_groupname_default)
        if len(default_collections) > 0:
            return default_collections[0]
        else:
            return None

    def has_default_metadata_collection(self, omsi_object):
        """
        Check whether the omsi API object (or h5py.Group) contains any a metadata collection with
        the default name.

        :param omsi_object: The omsi API object to be checked

        :returns: bool
        """
        default_meta_collection = self.get_default_metadata_collection(omsi_object=omsi_object)
        return default_meta_collection is not None


class omsi_file_metadata_collection(omsi_file_common):
    """
    Class for managing method specific data.

    **Use of super():**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common`.
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the __init__ function calls
    super(...).__init__(manager_group) with a single  parameter indicating the
    parent group.

    **Inherited Instance Variables**

    :ivar managed_group: The group that is managed by this object
    :ivar name: Name of the managed group

    """
    @classmethod
    def ___create___(cls,
                     parent_group,
                     group_name=None,
                     metadata=None,
                     type_attr_value="omsi_file_metadata_collection",
                     version_attr_value=omsi_format_metadata_collection.current_version,
                     flush_io=True):
        """
        Add a metadata management group to the given parent_group and populate it with
        the given metadata. NOTE: If the group already exists, then the existing group
        will be reused and only updated with the given metadata.

        :param parent_group: The h5py.Group where the method object should be created
        :type parent_group: h5py.Group
        :param group_name: The name of the group to be created or None if the default name
            for the general purpose metadata group should be used.
        :type group_name: str, None
        :param metadata: Additional metadata to be added to the collection after creation
        :type metadata: None, omsi.shared.metadata_data.metadata_value, omsi_shared.metadata_data.metadata_dict
        :param type_attr_value: Value to be used for the type attribute. This is used be derived classes
            to set their own type
        :type type_attr_value: str
        :param type_attr_value: Value to be used for the version attribute. This is used be derived classes
            to set their own version.
        :type version_attr_value: str
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 buffers are
                         flushed so that all data has been written to file

        :returns: `omsi_file_metadata_collection` object
        """
        metadata_group_name = group_name if \
            group_name is not None else \
            omsi_format_metadata_collection.metadata_collection_groupname_default
        # out_metadata_collection = None   # Created in both if and else block
        if metadata_group_name in parent_group.keys():
            out_metadata_collection = omsi_file_metadata_collection(metadata_group=parent_group[metadata_group_name])
            out_metadata_collection.add_metadata(metadata)
            if flush_io:
                parent_group.file.flush()
        else:
            metadata_group = parent_group.require_group(metadata_group_name)
            metadata_group.attrs[omsi_format_common.type_attribute] = type_attr_value
            metadata_group.attrs[omsi_format_common.version_attribute] = version_attr_value
            metadata_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
            out_metadata_collection = omsi_file_metadata_collection.__create_metadata_collection__(
                metadata_group=metadata_group,
                metadata=metadata)
        if flush_io:
            parent_group.file.flush()
        return out_metadata_collection

    @classmethod
    def __create_metadata_collection__(cls, metadata_group, metadata=None):
        """
        Add information about the metadata collection

        This function populates the group with the given metadata using
        add_metadata function and then creates and returns the instance
        of omsi_file_metadata_collection that manages the given group

        :returns: `omsi_file_metadata_collection` object responsible for managing the metadata_group.
        """
        new_metadata_collection = omsi_file_metadata_collection(metadata_group=metadata_group)
        new_metadata_collection.add_metadata(metadata=metadata)
        return new_metadata_collection

    def __init__(self, metadata_group):
        """
            Initialize the metadata collection object given the h5py object of the metadata collection

            :param metadata_group: The h5py object with the metadata collection group of the omsi hdf5 file.
            """
        super(omsi_file_metadata_collection, self).__init__(metadata_group)

    def __contains__(self, item):
        """
        Check if the given key is in the metadata
        :param item: Key to be checked
        :return: bool
        """
        return item in self.managed_group

    def keys(self):
        """
        Get a list of all metadata keys

        :return: List of string with the metadata keys
        """
        return self.managed_group.keys()

    def values(self):
        """
        Convenience function returning a list of all metadata objects.
        This is equivilant get_metadata(key=None).values(), however, for consistency with
        other dict-like interfaces this function returns a list of
        `omsi.shared.metadata_data.metadata_value` objects rather than
        the `omsi.shared.metadata_data.metadata_dict`

        :return: List of `omsi.shared.metadata_data.metadata_value` with all metadata
        """
        return self.get_metadata(key=None).values()

    def get_metadata(self, key=None):
        """
        Get dict with the full description of the metadata for the given key or all
        metadata if no key is given.

        :param key: The name of the metadata object to be retrieved. Default is None in
                    which case all metadata will be retrieved.

        :returns: `omsi.shared.metadata_data.metadata_value` object if a key is given
            or a `omsi.shared.metadata_data.metadata_dict` with all metadata
            if key is set to None.

        :raises: KeyError is raised in case that the specified key does not exist
        """
        descr_attr = omsi_format_metadata_collection.description_value_attribute
        unit_attr = omsi_format_metadata_collection.unit_value_attribute
        ontology_attr = omsi_format_metadata_collection.ontology_value_attribute
        isjson_attr = omsi_format_metadata_collection.is_json_dict_attribute
        if key is None:
            output_meta_dict = metadata_dict()
            for metadata_name, metadata_dataset in self.managed_group.iteritems():
                unit = None if unit_attr not in metadata_dataset.attrs else metadata_dataset.attrs[unit_attr]
                description = None if descr_attr not in metadata_dataset.attrs else metadata_dataset.attrs[descr_attr]
                ontology = None if ontology_attr not in metadata_dataset.attrs else \
                    json.loads(metadata_dataset.attrs[ontology_attr])
                md_value = metadata_dataset[:] if len(metadata_dataset.shape) > 0 else metadata_dataset[()]
                if isjson_attr in metadata_dataset.attrs:
                    try:
                        md_value = json.loads(md_value)
                    except:
                        log_helper.error(__name__, "Parsing of json metadata object failed for " + str(metadata_name))
                output_meta_dict[metadata_name] = metadata_value(
                    name=metadata_name,
                    value=md_value,
                    description=description,
                    unit=unit,
                    ontology=ontology)
            return output_meta_dict
        else:
            metadata_dataset = self.managed_group[key]
            unit = None if unit_attr not in metadata_dataset.attrs else metadata_dataset.attrs[unit_attr]
            description = None if descr_attr not in metadata_dataset.attrs else metadata_dataset.attrs[descr_attr]
            ontology = None if ontology_attr not in metadata_dataset.attrs else \
                json.loads(metadata_dataset.attrs[ontology_attr])
            return metadata_value(
                name=key,
                value=metadata_dataset[:],
                description=description,
                unit=unit,
                ontology=ontology)

    def add_metadata(self, metadata):
        """
        Add a new metadata entry

        :param metadata: Instance of `omsi.shared.metadata_data.metadata_value` or
            describing `omsi.shared.metadata_data.metadata_dict` with the
            metadata to be added.

        """
        if metadata is None:
            return

        if isinstance(metadata,  metadata_dict):
            for metadata_val in metadata.values():
                self.add_metadata(metadata_val)
            return

        metadata_dataset = None
        if isinstance(metadata['value'], basestring):
            metadata_dataset = self.managed_group.require_dataset(name=unicode(metadata['name']),
                                                                  shape=(1,),
                                                                  dtype=omsi_format_common.str_type)
            metadata_dataset[0] = unicode(metadata['value']) if \
                omsi_format_common.str_type_unicode else \
                str(metadata['value'])
        else:
            md_value = metadata['value']
            if md_value is not None:
                # Save different data types
                if isinstance(md_value, dict):
                    try:
                        json_value = json.dumps(md_value)
                        self.managed_group[metadata['name']] = json_value
                        isjson_attr_name = omsi_format_metadata_collection.is_json_dict_attribute
                        self.managed_group[metadata['name']].attrs[isjson_attr_name] = True
                        metadata_dataset = self.managed_group[metadata['name']]
                    except:
                        log_helper.error(__name__, "Metadata object for " + str(metadata['name'] +
                                         "is not JSON serializable"))
                        metadata_dataset = None

                elif isinstance(md_value, list):
                    try:
                        self.managed_group[metadata['name']] = np.asarray(metadata['value'])
                        metadata_dataset = self.managed_group[metadata['name']]
                    except:
                        log_helper.error(__name__, "Conversion and saving of list metadata failed for: " +
                                         str(metadata['name']))
                        metadata_dataset = None
                else:
                    self.managed_group[metadata['name']] = metadata['value']
                    metadata_dataset = self.managed_group[metadata['name']]
            else:
                self.managed_group[metadata['name']] = float_nan
                metadata_dataset = self.managed_group[metadata['name']]

        if metadata_dataset is not None:
            metadata_dataset.attrs[omsi_format_metadata_collection.description_value_attribute] = metadata['description']
            if metadata['unit'] is not None:
                metadata_dataset.attrs[omsi_format_metadata_collection.unit_value_attribute] = metadata['unit']
            if metadata['ontology'] is not None:
                ontology_value = metadata['ontology'] if isinstance(metadata['ontology'], basestring) \
                    else json.dumps(metadata['ontology'])
                metadata_dataset.attrs[omsi_format_metadata_collection.ontology_value_attribute] = ontology_value
