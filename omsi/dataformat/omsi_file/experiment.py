"""
OMSI file module for management of experiment data.
"""

import time

from omsi.dataformat.omsi_file.format import omsi_format_common, \
    omsi_format_experiment, \
    omsi_format_msidata
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager
from omsi.dataformat.omsi_file.msidata import omsi_file_msidata
from omsi.dataformat.omsi_file.methods import omsi_methods_manager
from omsi.dataformat.omsi_file.instrument import omsi_instrument_manager
from omsi.dataformat.omsi_file.analysis import omsi_analysis_manager
from omsi.dataformat.omsi_file.msidata import omsi_msidata_manager
import os


class omsi_experiment_manager(omsi_file_object_manager):
    """
    Experiment manager helper class used to define common functionality needed for
    experiment-related data. Usually, a class that defines a format that contains an
    omsi_file_experiment object will inherit from this class (in addition to omsi_file_common)
    to acquire the common features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar experiment_parent: The h5py.Group parent object containing the instrument object to be managed.

    """
    def __init__(self, experiment_parent):
        super(omsi_experiment_manager, self).__init__(experiment_parent)
        self.experiment_parent = experiment_parent

    def create_experiment(self, exp_identifier=None, flush_io=True):
        """
        Create a new group in the file for a new experiment and return the omsi_file_experiment
        object for the new experiment.

        :param exp_identifier: The string used to identify the analysis
        :type exp_identifier: string or None (default)
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that
                           all data has been written to file.

        :returns: omsi_file_experiment object for the newly created group for the experiment
        """
        return omsi_file_experiment.__create__(parent_group=self.experiment_parent,
                                               experiment_identifier=exp_identifier,
                                               experiment_index=None,
                                               flush_io=flush_io)

    @staticmethod
    def get_experiment_path(exp_index=None):
        """
        Based on the index of the experiment return the full path to the hdf5
        group containing the data for an experiment.

        :param exp_index: The index of the experiment.

        :returns: String indicating the path to the experiment.
        """

        expindexstr = ""
        if exp_index is not None:
            expindexstr = str(exp_index)
        return "/" + omsi_format_experiment.exp_groupname + expindexstr

    def get_num_experiments(self):
        """
        Get the number of experiments in this file.

        :returns: Integer indicating the number of experiments.
        """

        return omsi_file_common.get_num_items(self.experiment_parent, omsi_format_experiment.exp_groupname)

    def get_experiment(self, exp_index):
        """
        Get the omsi_format_experiment object for the experiment with the given index

        :param exp_index: The index of the requested experiment
        :type exp_index: uint

        :returns: h5py reference to the experiment with the given index. Returns None in case \
                  the experiment does not exist.
        """
        if exp_index < self.get_num_experiments():
            experiment_name = unicode(omsi_format_experiment.exp_groupname + str(exp_index))
            return omsi_file_experiment(self.experiment_parent[experiment_name])
        else:
            return None

    def get_experiment_by_identifier(self, exp_identifier_string):
        """
        Get the omsi_format_experiment object for the experiment with the given identifier.

        :param exp_identifier_string: The string used to identify the analysis
        :type exp_identifier_string: string

        :returns: Returns h5py object of the experiment group or None in case the experiment is not found.
        """
        # Iterate through all groups of the root folder
        for root_item in self.experiment_parent.items():
            if root_item[0].startswith(omsi_format_experiment.exp_groupname):
                cur_exp_id = omsi_file_experiment(self.experiment_parent[root_item[0]]).get_experiment_identifier()
                if cur_exp_id is not None:
                    if cur_exp_id[0] == exp_identifier_string:
                        return omsi_file_experiment(self.experiment_parent[root_item[0]])
        return None


class omsi_file_experiment(omsi_methods_manager,
                           omsi_instrument_manager,
                           omsi_analysis_manager,
                           omsi_msidata_manager,
                           omsi_file_common):
    """
    Class for managing experiment specific data

    **Use of super():**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common`.
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the __init__ function calls
    super(...).__init__(manager_group) with a single  parameter indicating the
    parent group.

    **Inherited instance variable:**

    :ivar managed_group: The group that is managed by this object
    :ivar methods_parent: The parent group containing the methods object (same as managed_group)
    :ivar instrument_parent: The parent group containing the instrument object (same as managed_group)
    :ivar name: Name of the managed group
    """

    @classmethod
    def __create__(cls,
                   parent_group,
                   experiment_identifier,
                   experiment_index=None,
                   flush_io=True):
        """
        Create an an experiment group. Populate the group using the omsi_file_experiment.__populate_experiment__
        function.

        :param parent_group: The parent h5py group where the experiment group should be created in.
        :type parent_group. h5py.Group
        :param experiment_identifier: Identifier string for the experiment
        :type experiment_identifier: string, None
        :param experiment_index: Specify the index of the experiment. Use None to generate the next experiment.
        :type experiment_index: uint or None
        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed so that all
                         data has been written to file

        :returns: The function returns the omsi_file_experiment object responsible for managing the group

        :raises: IndexError is raised in case that the experiment with the given index already exists.

        """
        if experiment_index is None:
            experiment_index = omsi_file_common.get_num_items(parent_group,
                                                              omsi_format_experiment.exp_groupname)
        experiment_name = omsi_format_experiment.exp_groupname + str(experiment_index)
        # Using require_group ensures that the group is not overwritten if it
        # already exists
        if experiment_name in parent_group.keys():
            raise IndexError("The experiment with the given index already exists")
        experiment_group = parent_group.require_group(experiment_name)
        experiment_group.attrs[omsi_format_common.type_attribute] = "omsi_file_experiment"
        experiment_group.attrs[omsi_format_common.version_attribute] = omsi_format_experiment.current_version
        experiment_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())

        # Populate the group and create the omsi_file_experiment object
        experiment_object = omsi_file_experiment.__populate_experiment__(experiment_group=experiment_group,
                                                                         experiment_identifier=experiment_identifier)
        # Flush I/O if needed
        if flush_io:
            parent_group.file.flush()
        # Return the omsi_file_experiment object
        return experiment_object

    @classmethod
    def __populate_experiment__(cls,
                                experiment_group,
                                experiment_identifier):
        """
        Populate the given experiment group.

        :param experiment_group: The group to be populated.
        :param experiment_identifier: The experiment identifier to be used.
        :type experiment_identifier: string or None. If None then 'undefined' will be used.
        :return:
        """
        experiment_identifier_dataset = experiment_group.require_dataset(name=unicode(
            omsi_format_experiment.exp_identifier_name), shape=(1,), dtype=omsi_format_common.str_type)
        if experiment_identifier is not None:
            if omsi_format_common.str_type_unicode:
                experiment_identifier_dataset[0] = experiment_identifier
            else:
                experiment_identifier_dataset[0] = str(experiment_identifier)
        else:
            experiment_identifier_dataset[0] = "undefined"
        return omsi_file_experiment(experiment_group)

    def __init__(self, exp_group):
        """
        Initalize the experiment object given the h5py object of the experiment group

        :param exp_group: The h5py object with the experiment group of the omsi hdf5 file.
        """
        super(omsi_file_experiment, self).__init__(exp_group)
        # The following initialization are performed by the super call
        # self.managed_group = exp_group
        # self.methods_parent = self.managed_group
        # self.instrument_parent = self.managed_group
        # self.name = self.managed_group.name

    def get_experiment_index(self):
        """
        Determine the index of the experiment based on the name of the group
        :return: Integer index of the experiment
        """
        return int(os.path.basename(self.managed_group.name).lstrip(omsi_format_experiment.exp_groupname))

    ##########################################################
    #  Access the data directly associated with the experiment
    ##########################################################
    def get_experiment_identifier(self):
        """
        Get the HDF5 dataset with the identifier description for the experiment.

        :returns: h5py object of the experiment identifier or None in case not present
        """
        try:
            return self.managed_group[unicode(omsi_format_experiment.exp_identifier_name)]
        except KeyError:
            return None

    def set_experiment_identifier(self, identifier):
        """
        Overwrite the current identfier string for the experiment with the given string

        :param identifier: The new experiment identifier string.
        """
        # Get the name of the intrument
        expid = self.get_experiment_identifier()
        # Create the dataset for the id name if it does not exist
        if expid is None:
            expid = self.managed_group.require_dataset(name=unicode(
                omsi_format_experiment.exp_identifier_name), shape=(1,), dtype=omsi_format_common.str_type)

        if omsi_format_common.str_type_unicode:
            expid[0] = identifier
        else:
            expid[0] = str(identifier)

    ###########################################################
    # Get sub-group object associated with the experiment
    ###########################################################
    def get_instrument_info(self, check_parent=False):
        """
        Inherited from omsi_instrument_manager parent class. Overwritten here to change the default
        parameter setting for check_parent. See `omsi.dataformat.omsi_file.instrument.omsi_instrument_manager`
        for details.
        """
        return super(omsi_file_experiment, self).get_instrument_info(check_parent=check_parent)

    def get_method_info(self, check_parent=False):
        """
        Inherited from omsi_method_manager parent class. Overwritten here to change the default
        parameter setting for check_parent. See `omsi.dataformat.omsi_file.methods.omsi_method_manager`
        for details
        """
        return super(omsi_file_experiment, self).get_method_info(check_parent=check_parent)
