"""
Base module for managing of dependencies between data in OpenMSI HDF5 files
"""

import time

from omsi.dataformat.omsi_file.format import omsi_format_common, \
    omsi_format_dependencies, \
    omsi_format_dependencydata
from omsi.dataformat.omsi_file.common import omsi_file_common, omsi_file_object_manager


class omsi_dependencies_manager(omsi_file_object_manager):
    """
    Dependencies manager helper class
    used to define common functionality needed for managing dependencies.
    Usually, a class that defines a format that contains an omsi_file_dependencies object
    will inherit from this class (in addition to omsi_file_common) to acquire the common
    features.

    For more details see: `omsi.dataforamt.omsi_file.omsi_common.omsi_file_object_manager`

    :ivar dependencies_parent: h5py.Group object containing the dependencies object(s) to be managed
    :ivar dependencies: omsi_file_dependencies object managed by this object or None

    """
    def __init__(self, dependencies_parent):
        """

        :param dependencies_parent: Parent group containing the dependencies object to be managed
        """
        super(omsi_dependencies_manager, self).__init__(dependencies_parent)
        self.dependencies_parent = dependencies_parent
        try:
            self.dependencies = omsi_file_dependencies(
                self.dependencies_parent[unicode(omsi_format_dependencies.dependencies_groupname)])
        except:
            self.dependencies = None

    def create_dependencies(self, dependencies_data_list=None):
        """
        Create a managed group for storing data dependencies

        :param parent_group: The h5py.Group for which the dependencies group should be initialized.
        :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies.
               Default is None which is mapped to an empty list []

        :returns: omsi_file_dependencies object created by the function.
        """
        return omsi_file_dependencies.__create__(parent_group=self.dependencies_parent,
                                                 dependencies_data_list=dependencies_data_list)

    def add_dependency(self, dependency, flush_io=True):
        """
        Create a new dependency for this dataset

        :param dependency: omsi.shared.omsi_dependency object describing the data dependency

        :param flush_io: Call flush on the HDF5 file to ensure all HDF5 bufferes are flushed
                         so that all data has been written to file

        :returns: omsi_file_dependencydata object with the dependency data or None in case that
                  an error occurred and the dependency has not been generated.
        """
        if dependency is None:
            raise ValueError('dependency data object may not be None.')
        if not self.dependencies:
            self.dependencies = omsi_file_dependencies.__create__(parent_group=self.dependencies_parent,
                                                                  dependencies_data_list=None)
        out_dependency = self.dependencies.add_dependency(dependency)
        if flush_io:
            self.dependencies_parent.file.flush()
        return out_dependency

    def get_all_dependency_data(self, omsi_dependency_format=True):
        """
        Get all direct dependencies associated with the data object.

        :param omsi_dependency_format: Should the dependencies be retrieved as omsi_analysis_dependency object (True)
                                       or as an omsi_file_dependencydata object (False).

        :returns: List omsi_dependency objects containing either omsi file API objects or h5py objects for the
                  dependencies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def get_all_dependency_data_recursive(self, omsi_dependency_format=True):
        """
        Get all direct and indirect dependencies associated with the data object.

        NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

        :param omsi_dependency_format: Should the dependencies be retrieved as omsi_dependency object (True)
               or as an omsi_file_dependencydata object (False)

        :returns: List omsi_analysis_data objects containing either omsi file API interface objects or h5py
                  objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data_recursive(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def has_dependencies(self):
        """
        Check whether any dependencies exists for this datasets
        """
        if self.dependencies:
            return True
        else:
            return False


class omsi_file_dependencies(omsi_file_common):
    """
    Class for managing collections of dependencies.

    ** Use of super()**

    This class inherits from `omsi.dataformat.omsi_file.common.omsi_file_common`.
    Consistent with the design pattern for multiple inheritance of the
    `omsi.dataformat.omsi_file` module, the __init__ function calls
    super(...).__init__(manager_group) with a single  parameter indicating the
    parent group.

    """
    @classmethod
    def __create__(cls, parent_group, dependencies_data_list=None):
        """
        Create a managed group for storing data dependencies

        :param parent_group: The h5py.Group for which the dependencies group should be initialized.
        :param dependencies_data_list: List of omsi_dependency objects to be stored as dependencies.
               Default is None which is mapped to an empty list []

        :returns: omsi_file_dependencies object created by the function.
        """
        dependencies_group = parent_group.require_group(omsi_format_dependencies.dependencies_groupname)
        dependencies_group.attrs[omsi_format_common.type_attribute] = "omsi_file_dependencies"
        dependencies_group.attrs[omsi_format_common.version_attribute] = omsi_format_dependencies.current_version
        dependencies_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        dependencies_object = omsi_file_dependencies.__populate_dependencies__(
            dependencies_group=dependencies_group,
            dependencies_data_list=dependencies_data_list)
        return dependencies_object

    @classmethod
    def __populate_dependencies__(cls, dependencies_group, dependencies_data_list=None):
        """
        Populate the given h5py group with the dependencies data.

        NOTE: This is a private helper function. Use the corresponding create_dependencies function \
        of omsi_file_analysis and omsi_file_msidata to create a completely new collection of dependencies.

        :param dependencies_group: h5py.Group object that should be used to store the dependencies data.
        :type dependencies_group: h5py.Group
        :param dependencies_data_list: List of dependencies to be stored
        :type dependencies_data_list: List of omsi.shared.omsi_dependency_data objects to be stored.
        """
        omsiobj = omsi_file_dependencies(dependencies_group)
        if dependencies_data_list:
            for dependency_obj in dependencies_data_list:
                omsiobj.add_dependency(dependency_data=dependency_obj)
        return omsiobj

    def __init__(self, dependencies_group):
        super(omsi_file_dependencies, self).__init__(dependencies_group)
        self.managed_group = dependencies_group
        self.name = self.managed_group.name

    def add_dependency(self, dependency_data):
        """
        Add a new dependency to the collection.

        :param dependency_data: The analysis dependency specification.
        :type dependency_data: omsi.shared.omsi_dependency_data

        :returns: the newly created omsi_file_dependencydata object
        :raises: KeyError in case that a dependency with the same name already exists
        """
        return omsi_file_dependencydata.__create__(parent_group=self.managed_group,
                                                   dependency_data=dependency_data)

    def get_omsi_file_dependencydata(self, name):
        """
        Retrieve the omsi_file_dependencydata object for the dependency with the given name.
        """
        return omsi_file_common.get_omsi_object(self.managed_group[name])

    def get_dependency_omsiobject(self, name, recursive=True):
        """
        Get the omsi file API object corresponding to the object the dependency is pointing to.

       :param name: Name of the dependency opbject to be loaded .
       :param recursive: Should the dependency be resolved recursively, i.e., if the dependeny points
                to another dependencies. Default=True.

       :returns: An omsi file API object (e.g., omsi_file_analysis or omsi_file_msidata) if the link points
                to a group or the h5py.Dataset the link is pointing to.
        """
        return self.get_omsi_file_dependencydata(name).get_dependency_omsiobject(recursive=recursive)

    def get_all_dependency_data(self, omsi_dependency_format=True):
        """
        Get all direct dependencies associdated with the analysis.

        :param omsi_dependency_format: Should the dependcies be retrieved as omsi_analysis_dependency object
                (True) or as an omsi_file_dependencydata object (False).

        :returns: List omsi_dependency objects containing either omsi file API objects or h5py objects for
                the dependcies. Access using [index]['name'] and [index]['data'].
        """
        output_list = []
        for item_obj in self.managed_group.items():
            omsi_object = omsi_file_dependencydata(self.managed_group[unicode(item_obj[0])])
            if omsi_dependency_format:
                output_list.append(omsi_object.get_omsi_dependency())
            else:
                output_list.append(omsi_object)
        return output_list

    def get_all_dependency_data_recursive(self, omsi_dependency_format=True):
        """
        Get all direct and indirect dependencies associdated with the analysis.

        NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

        :param omsi_dependency_format: Should the dependencies be retrieved as omsi_dependency object (True)
                or as an omsi_file_dependencydata object (False)

        :returns: List omsi_analysis_data objects containing either omsi file API interface objects or
                h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        output_list = []
        for item_obj in self.managed_group.items():
            omsi_obj = omsi_file_common.get_omsi_object(self.managed_group[item_obj[0]])
            # If we can have recusive dependencies
            if isinstance(omsi_obj, omsi_file_analysis):
                it_depend = omsi_obj.get_all_dependency_data_recursive(
                    omsi_dependency_format=omsi_dependency_format)
                output_list = output_list + it_depend
            try:
                omsi_object = omsi_file_dependencydata(
                    self.managed_group[unicode(item_obj[0])])
                if omsi_dependency_format:
                    output_list.append(omsi_object.get_omsi_dependency())
                else:
                    output_list.append(omsi_object)
            except:
                import sys
                print "WARNING: Error occured in omsi_file_dependencies::get_all_dependency_data_recursive(...):  " + \
                      unicode(item_obj[0]) + "   :" + str(sys.exc_info())
        return output_list

    def get_all_dependency_data_graph(self,
                                      include_omsi_dependency=False,
                                      include_omsi_file_dependencydata=False,
                                      recursive=True,
                                      level=0,
                                      name_key='name',
                                      prev_nodes=None,
                                      prev_links=None,
                                      parent_index=None,
                                      metadata_generator=None,
                                      metadata_generator_kwargs=None):
        """
        Get all direct and indirect dependencies associdated with the analysis in form of a graph describing
        all nodes and links in the provenance hierarchy.

        NOTE: Circular dependcies are not supported and should not be possible in practice anyways.

        :param include_omsi_dependency: Should the omsi_dependency object be included in the entries
                    in the nodes dict?
        :param include_omsi_file_dependencydata: Should the omsi_file_dependencydata object be included in
                    the entries in the nodes dict?
        :param recursive: Should we trace dependencies recursively to construct the full graph, or only the
                    direct dependencies. Default true (ie., trace recursively)
        :param name_key: Which key should be used in the dicts to indicate the name of the object?
                    Default value is 'name'
        :param level: Integer used to indicated the recursion level. Default value is 0.
        :param prev_nodes: List of nodes that have been previously generated. Note, this list will be
                    modified by the call. Note, each node is represented by a dict which is expected to
                    have at least the following keys defined, path, name_key, level (name_key refers to the key
                    defined by the input parameter name_key.
        :param prev_links: Previouly established links in the list of nodes. Note, this list will be
                    modified by the call.
        :param parent_index: Index of the parent node in the list of prev_nodes for this call.
        :param metadata_generator: Optional parameter. Pass in a function that generates additional metadata
                    about a given omsi API object. Note, the key's level and path and name (i.e., name_key) are
                    already set by this function. The metadata_generator may overwrite these key's, however,
                    the path has to be unique as it is used to identify duplicate nodes. Overwriting the path
                    with a non-unique value, hence, will lead to errors (missing entries) when generating the graph.
                    Note, the metadata_generator function needs to support the following keyword arguments:

                          * inDict : The dictionary to which the metadata should be added to.
                          * obj : The omsi file API object for which metadata should be generated
                          * name : A qualifying name for the object
                          * name_key : The key to be used for storing the name
        :param metadata_generator_kwargs: Dictionary of additional keyword arguments that should be passed to
                    the metadata_generator function.

        :returns: Dictionary containing two lists. 1) nodes : List of dictionaries, describing the elements
                  in the dependency graph. 2) links : List of tuples with the links in the graph. Each
                  tuple consists of two integer indices for the nodes list. For each node the following
                  entries are given:

                * omsi_dependency: Optional key used to store the corresponding omsi_dependency object.
                        Used only of include_omsi_dependency is True.
                * omsi_file_dependencydata: Optional key used to store the corresponding
                        omsi_file_dependencydata object. Used only of include_omsi_file_dependencydata is True.
                * name : Name of the dependency. The actual key is sepecified by name_key
                * level : The recursion level at which the object occurs.
                * ... : Any other key/value pairs from the omsi_dependency dict.

        """
        from omsi.dataformat.omsi_file.analysis import omsi_file_analysis
        from omsi.dataformat.omsi_file.msidata import omsi_file_msidata

        import os
        if prev_nodes:
            nodes = prev_nodes
        else:
            nodes = []
        if prev_links:
            links = prev_links
        else:
            links = []
        if metadata_generator_kwargs:
            metagenkwargs = metadata_generator_kwargs
        else:
            metagenkwargs = {}

        def find_node_index(path_string):
            """
            Internal helper funtion used to check whether a node already exists in the graph

            :param path_string: The path to object used for comparison.
            """
            if prev_nodes:
                for i in range(len(prev_nodes)):
                    if nodes[i]['path'] == path_string:
                        return i
            return None

        # Iterate through all dependencies
        for item_obj in self.managed_group.items():
            try:
                # 1) Check if a node exists that represents the current
                # dependency
                dep_obj = omsi_file_dependencydata(
                    self.managed_group[unicode(item_obj[0])])
                omsi_obj = dep_obj.get_dependency_omsiobject(recursive=True)
                curr_path = omsi_obj.name
                curr_name = os.path.basename(curr_path)
                curr_index = find_node_index(omsi_obj.name)
                # 2) If no node exists for the current dependency then add a new node
                if not curr_index:
                    # 2.1)Add a new node with the basic metadata
                    nodes.append({})
                    curr_index = len(nodes) - 1
                    curr_node = nodes[curr_index]
                    curr_node['level'] = level
                    curr_node[name_key] = curr_name
                    curr_node['path'] = curr_path
                    # print curr_path
                    if include_omsi_dependency:
                        curr_node[
                            'omsi_dependency'] = dep_obj.get_omsi_dependency()
                    if include_omsi_file_dependencydata:
                        curr_node['omsi_file_dependencydata'] = dep_obj
                    # 2.2) Expand the metadata dict for the current object with
                    # any user-defined metadata
                    if metadata_generator:
                        metadata_generator(inDict=curr_node,
                                           obj=omsi_obj,
                                           name=os.path.basename(omsi_obj.name),
                                           name_key=name_key,
                                           **metagenkwargs)

                # print "PARENT:"+str(parent_index)
                # 3) Add a link from the current node to its parent
                if parent_index is not None:
                    templink = {'source': parent_index, 'target': curr_index}
                    links.append(templink)
                    # print "LINK: "+str(templink)

                # 4) If we can have recusive dependencies that we should trace,
                # then trace them back and add them to the nodes and links
                # lists
                if recursive:
                    if isinstance(omsi_obj, omsi_file_analysis) or isinstance(omsi_obj, omsi_file_msidata):
                        if omsi_obj.dependencies:
                            omsi_obj.dependencies.get_all_dependency_data_graph(
                                include_omsi_dependency=include_omsi_dependency,
                                include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                                level=level +
                                1,
                                name_key=name_key,
                                prev_nodes=nodes,
                                prev_links=prev_links,
                                parent_index=curr_index,
                                metadata_generator=metadata_generator,
                                metadata_generator_kwargs=metagenkwargs
                            )
                # print "-----B----"+str(links)
            except:
                import sys
                print "WARNING: Error occured in omsi_file_dependencies::get_all_dependency_data_graph(...):" \
                      + str(sys.exc_info())

        return nodes, links


##########################################################################
##########################################################################
#   The main class for accessing dependency datasets                           #
##########################################################################
##########################################################################
class omsi_file_dependencydata(omsi_file_common):

    """Class for managing data groups used for storing data dependencies"""
    @classmethod
    def __create__(cls, parent_group, dependency_data):
        """
        Create a new dependency group within the given parent group and populate it
        with the required datasets using omsi_file_dependencydata.__populate_dependency__
        based on the given dependency information.

        :param dependency_group: The h5py group opbject to which the dependency dataset should be added
        :type dependency_group: h5py.Group
        :param dependency_data: The analysis dependency specification.
        :type dependency_data: omsi.shared.omsi_dependency_data

        :returns: omsi_file_dependencydata object for management of the newly created object

        """
        dependency_name = dependency_data['link_name']
        # Check if the group already exists and raise an error if it does
        if dependency_name in parent_group.keys():
            raise KeyError("Another dependencies with the same index already exists.")
        new_dep_group = parent_group.require_group(dependency_name)
        new_dep_group.attrs[omsi_format_common.type_attribute] = "omsi_file_dependencydata"
        new_dep_group.attrs[omsi_format_common.version_attribute] = omsi_format_dependencydata.current_version
        new_dep_group.attrs[omsi_format_common.timestamp_attribute] = str(time.ctime())
        dependency_object = omsi_file_dependencydata.__populate_dependency__(dependency_group=new_dep_group,
                                                                             dependency_data=dependency_data)
        return dependency_object

    @classmethod
    def __populate_dependency__(cls, dependency_group, dependency_data):
        """
        Populate the given dependency_group with the associated data.

        :param dependency_group: The h5py group opbject to which the dependency dataset should be added
        :type dependency_group: h5py.Group
        :param dependency_data: The analysis dependency specification.
        :type dependency_data: omsi.shared.omsi_dependency

        :returns: omsi_file_dependencydata object for management of the newly created object

        """
        dep_group = dependency_group
        # Save the name of the parameter
        param_name_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_parameter), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            param_name_data[0] = dependency_data['param_name']
        else:
            param_name_data[0] = str(dependency_data['param_name'])
        # Save the selection
        selection_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_selection), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['selection'] is not None:
            from omsi.shared.omsi_data_selection import check_selection_string
            # This should always be True since omsi_dependency checks for this
            # but we need to be sure.
            if check_selection_string(dependency_data['selection']):
                if omsi_format_common.str_type_unicode:
                    selection_data[0] = dependency_data['selection']
                else:
                    selection_data[0] = str(dependency_data['selection'])
            else:
                selection_data[0] = ""
                errormessage = "Invalid selection string given for data dependency : " + \
                    str(dependency_data['selection'])
                raise ValueError(errormessage)
        else:
            selection_data[0] = ""
        # Save the main omsi object
        mainname_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_mainname), shape=(1,), dtype=omsi_format_common.str_type)
        if omsi_format_common.str_type_unicode:
            mainname_data[0] = unicode(dependency_data['omsi_object'].name)
        else:
            mainname_data[0] = str(dependency_data['omsi_object'].name)
        # Save the additional dataset name
        dataset_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_datasetname), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['dataname']:
            if omsi_format_common.str_type_unicode:
                dataset_data[0] = dependency_data['dataname']
            else:
                dataset_data[0] = str(dependency_data['dataname'])
        else:
            dataset_data[0] = u''

        return omsi_file_dependencydata(dep_group)

    def __init__(self, dependency_group):
        """
        Create a new omsi_file_dependencydata object for the given h5py.Group

        :param dependency_group: h5py.Group object with the dependency data
        :type dependency_group: h5py.Group with a corresponding omsi type
        """
        super(omsi_file_dependencydata, self).__init__(dependency_group)
        self.managed_group = dependency_group
        self.name = self.managed_group.name

    def __getitem__(self, key):
        """
        Retrieve the h5py dataset for any of the dependency datasetes or the omsi file API object.
        If the key is a numpy selection (slicing) then the data is retrieved from the link the
        dependency is pointing so that we can interact with the dependency as if it were the
        linked object.
        """
        if isinstance(key, str) or isinstance(key, unicode):
            if key == omsi_format_dependencydata.dependency_mainname or \
               key == omsi_format_dependencydata.dependency_datasetname:
                return omsi_file_common.get_omsi_object(self.managed_group[key])
            else:
                return self.managed_group[key][0]
        else:
            omsi_object_name = self.managed_group[
                unicode(omsi_format_dependencydata.dependency_mainname)][0]
            h5py_object = self.managed_group.file[unicode(omsi_object_name)]
            omsi_object = omsi_file_common.get_omsi_object(h5py_object)
            dataset_name = self.managed_group[
                unicode(omsi_format_dependencydata.dependency_datasetname)][0]
            if len(dataset_name) > 0:
                return omsi_object[dataset_name][key]
            else:
                return omsi_object[key]

    def get_dependency_type(self, recursive=True):
        """
        Indicated the type of the object the dependency is pointing to.

        :param recursive: Should dependencies be resolved recursively, i.e., if the dependency points to another
                         dependencies. Default=True.

        :returns: String indicating the class of the omsi file API class that is suited to manage the dependency
                  link or the name of the corresponding h5py class.
        """
        omsi_object = self.get_dependency_omsiobject(recursive=recursive)
        return omsi_object.__class__.__name__

    def get_selection_string(self):
        """
        String indicating the applied selection. This is an empty string in case no selection was applied.

        :returns: Selection string. See the omsi.shared.omsi_data_selection for helper functions to
                  deal with selection strings.
        """
        return self.managed_group[unicode(omsi_format_dependencydata.dependency_selection)][0]

    def get_parameter_name(self):
        """
        Get the string indicating the name of the dependend parameter of the analysis.

        :returns: String of the parameter name that has the dependency.
        """
        return self.managed_group[unicode(omsi_format_dependencydata.dependency_parameter)][0]

    def get_dataset_name(self):
        """
        Get the string indicating the name of dataset. This may be empty as it is only used if the
        dependency points to an objec within a managed omsi API object.

        :returns: String indicating the name of the optional dataset.
        """
        try:
            return self.managed_group[unicode(omsi_format_dependencydata.dependency_datasetname)][0]
        except:
            return ""

    def get_dependency_omsiobject(self, recursive=True):
        """
        Get the omsi file API object corresponding to the object the dependency is pointing to.

        :param recursive: Should dependencies be resolved recursively, i.e., if the dependency points
                          to another dependencies. Default=True.

        :returns: An omsi file API object (e.g., omsi_file_analysis or omsi_file_msidata) if the link
                  points to a group or the h5py.Dataset the link is pointing to.
        """
        omsi_object_name = self.managed_group[
            unicode(omsi_format_dependencydata.dependency_mainname)][0]
        h5py_object = self.managed_group.file[unicode(omsi_object_name)]
        omsi_object = omsi_file_common.get_omsi_object(h5py_object,
                                                       resolve_dependencies=recursive)
        return omsi_object

    def get_link_name(self):
        """
        Get the name of the dependency link

        :returns: String indicating the name of the dependency link.
        """
        return self.managed_group.name.split("/")[-1]

    def get_omsi_dependency(self):
        """
        Get the dependency information as an omsi.shared.omsi_dependency object
        (as defined in the omsi.shared.omsi_dependency module)

        :returns: omsi_dependency object with all the dependency data.
        """
        from omsi.shared.omsi_dependency import omsi_dependency
        output_dependency = omsi_dependency()
        try:
            output_dependency['param_name'] = self.get_parameter_name()
        except:
            output_dependency['param_name'] = None
        try:
            output_dependency['link_name'] = self.get_link_name()
        except:
            output_dependency['link_name'] = None
        try:
            output_dependency['selection'] = self.get_selection_string()
        except:
            output_dependency['selection'] = None
        try:
            output_dependency['omsi_object'] = self.get_dependency_omsiobject()
        except:
            output_dependency['omsi_object'] = None
        try:
            output_dependency['dataname'] = self.get_dataset_name()
        except:
            output_dependency['dataname'] = None

        return output_dependency
