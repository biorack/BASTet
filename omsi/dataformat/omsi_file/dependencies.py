"""
Base module for managing of dependencies between data in OpenMSI HDF5 files
"""

import time
import os
import warnings

import h5py

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
        Create a managed group for storing data dependencies if none exists and store
        the given set of dependencies in it. If a self.dependencies object already
        exists, then the given dependencies will be added.

        This is effectively a shortcut to omsi_file_dependencies.__create___(...) with specific
        settings for the current dependencies object managed by self.

        :param dependencies_data_list: List of dependency_dict objects to be stored as dependencies.
               Default is None which is mapped to an empty list []

        :returns: omsi_file_dependencies object created by the function.
        """
        # If we don't have any dependency group then create it and all dependencies
        if not self.dependencies:
            self.dependencies = omsi_file_dependencies.__create__(parent_group=self.dependencies_parent,
                                                                  dependencies_data_list=dependencies_data_list)
        # Otherwise add all dependencies to the current dependency group
        else:
            for dependency in dependencies_data_list:
                self.add_dependency(dependency=dependency)
        return self.dependencies

    def add_dependency(self, dependency, flush_io=True):
        """
        Create a new dependency for this dataset

        :param dependency: omsi.shared.dependency_dict object describing the data dependency

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

        This is convenience function providing access to
        self.dependencies.get_all_dependency_data(...) which is a function
        of omsi_file_dependencies class.

        :param omsi_dependency_format: Should the dependencies be retrieved as omsi_analysis_dependency object (True)
                                       or as an omsi_file_dependencydata object (False).

        :returns: List dependency_dict objects containing either omsi file API objects or h5py objects for the
                  dependencies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            return self.dependencies.get_all_dependency_data(omsi_dependency_format=omsi_dependency_format)
        else:
            return []

    def get_all_dependency_data_recursive(self,
                                          omsi_dependency_format=True,
                                          omsi_main_parent=None,
                                          dependency_list=None):
        """
        Get all direct and indirect dependencies associated with the data object.

        This is convenience function providing access to
        self.dependencies.get_all_dependency_data_recursive(...) which is a function
        of omsi_file_dependencies class.

        **NOTE:** omsi_main_parent and omsi_main_parent are used primarily to ensure that the
        case of circular dependencies are supported properly. Circular dependencies may
        occur in the case of semantic dependencies (rather than pure use dependencies), e.g.,
        two datasets that are related modalities may reference each other, e.g., MS1 pointing
        to related MS2 data and the MS2 datasets referencing the corresponding MS1 datasets.

        :param omsi_dependency_format: Should the dependencies be retrieved as dependency_dict object (True)
               or as an omsi_file_dependencydata object (False)

        :param omsi_main_parent: The main parent for which the dependencies are calculated. This is needed
            to avoid recursion back into the main parent for which we are computing dependencies and avoiding
            that it is added itself as a dependency for itself. If set to None, then we will use our own
            self.dependencies_parent object

        :param dependency_list: List of previously visited/created dependencies. This is needed only
            to avoid deep recursion and duplication due to circular dependencies

        :returns: List analysis_data objects containing either omsi file API interface objects or h5py
                  objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if self.dependencies is not None:
            if omsi_main_parent is None:
                omsi_main_parent = omsi_file_common.get_omsi_object(self.dependencies_parent)
            return self.dependencies.get_all_dependency_data_recursive(
                omsi_dependency_format=omsi_dependency_format,
                omsi_main_parent=omsi_main_parent,
                dependency_list=dependency_list)
        else:
            return []

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
        Get all direct and indirect dependencies associated with the analysis in form of a graph describing
        all nodes and links in the provenance hierarchy.

        This is convenience function providing access to
        self.dependencies.get_all_dependency_data_graph(...) which is a function
        of omsi_file_dependencies class.

        :param include_omsi_dependency: Should the dependency_dict object be included in the entries
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

                * dependency_dict: Optional key used to store the corresponding dependency_dict object.
                        Used only of include_omsi_dependency is True.
                * omsi_file_dependencydata: Optional key used to store the corresponding
                        omsi_file_dependencydata object. Used only of include_omsi_file_dependencydata is True.
                * name : Name of the dependency. The actual key is sepecified by name_key
                * level : The recursion level at which the object occurs.
                * ... : Any other key/value pairs from the dependency_dict dict.

        """
        if self.dependencies:
            return self.dependencies.get_all_dependency_data_graph(
                include_omsi_dependency=include_omsi_dependency,
                include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                recursive=recursive,
                level=level,
                name_key=name_key,
                prev_nodes=prev_nodes,
                prev_links=prev_links,
                parent_index=parent_index,
                metadata_generator=metadata_generator,
                metadata_generator_kwargs=metadata_generator_kwargs)
        else:
            empty_graph = {'nodes': [], 'links': []}
            return empty_graph

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
        :param dependencies_data_list: List of dependency_dict objects to be stored as dependencies.
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
        Get all direct dependencies associated with the analysis.

        :param omsi_dependency_format: Should the dependencies be retrieved as omsi_analysis_dependency object
                (True) or as an omsi_file_dependencydata object (False).

        :returns: List dependency_dict objects containing either omsi file API objects or h5py objects for
                the dependencies. Access using [index]['name'] and [index]['data'].
        """
        output_list = []
        for item_obj in self.managed_group.items():
            omsi_object = omsi_file_dependencydata(self.managed_group[unicode(item_obj[0])])
            if omsi_dependency_format:
                output_list.append(omsi_object.get_omsi_dependency())
            else:
                output_list.append(omsi_object)
        return output_list

    def get_all_dependency_data_recursive(self,
                                          omsi_dependency_format=True,
                                          omsi_main_parent=None,
                                          dependency_list=None):
        """
        Get all direct and indirect dependencies associated with the analysis.

        **NOTE:** omsi_main_parent and omsi_main_parent are used primarily to ensure that the
        case of circular dependencies are supported properly. Circular dependencies may
        occur in the case of semantic dependencies (rather than pure use dependencies), e.g.,
        two datasets that are related modalities may reference each other, e.g., MS1 pointing
        to related MS2 data and the MS2 datasets referencing the corresponding MS1 datasets.

        :param omsi_dependency_format: Should the dependencies be retrieved as dependency_dict object (True)
                or as an omsi_file_dependencydata object (False)

        :param omsi_main_parent: The main parent for which the dependencies are calculated. This is needed
            to avoid recursion back into the main parent for which we are computing dependencies and avoiding
            that it is added itself as a dependency for itself. If set to None, then we will use the omsi_object
            associated with the parent group of the dependency group.

        :param dependency_list: List of previously visited/created dependencies. This is needed only
            to avoid deep recursion and duplication due to circular dependencies

        :returns: List analysis_data objects containing either omsi file API interface objects or
                h5py objects for the dependcies. Access using [index]['name'] and [index]['data'].
        """
        if omsi_main_parent is None:
            omsi_main_parent = omsi_file_common.get_omsi_object(self.managed_group.parent)
        if dependency_list is None:
            dependency_list = []
        output_list = []
        for item_obj in self.managed_group.items():
            # omsi_obj is a omsi_file_dependencydata object
            omsi_obj = omsi_file_common.get_omsi_object(self.managed_group[item_obj[0]])
            # dependency_omsi_obj is the instance of the omsi API object we are depending on
            dependency_omsi_obj = omsi_obj.get_dependency_omsiobject()

            # Add the object to the list of dependencies
            try:
                omsi_object = omsi_file_dependencydata(self.managed_group[unicode(item_obj[0])])
                if omsi_dependency_format:
                    output_list.append(omsi_object.get_omsi_dependency())
                else:
                    output_list.append(omsi_object)
            except:
                import sys
                warnings.warn("WARNING: Error occurred in omsi_file_dependencies::get_all_dependency_data_recursive(...):  " + \
                      unicode(item_obj[0]) + "   :" + str(sys.exc_info()))

            # If we can have recursive dependencies then follow them
            if isinstance(dependency_omsi_obj, omsi_dependencies_manager):
                # Check if the same omsi_object is already in the list of dependencies. If yes,
                # then we do not need to do the recursion any more
                check_dependencies_recursively = dependency_omsi_obj != omsi_main_parent
                curr_full_dependency_list = dependency_list+output_list
                if omsi_dependency_format:
                    for prior_dependency in curr_full_dependency_list:
                        if prior_dependency['omsi_object'] == dependency_omsi_obj:
                            check_dependencies_recursively = False
                            break
                else:
                    for prior_dependency in curr_full_dependency_list:
                        if prior_dependency.get_dependency_omsiobject() == dependency_omsi_obj:
                            check_dependencies_recursively = False
                            break
                if check_dependencies_recursively:
                    it_depend = dependency_omsi_obj.get_all_dependency_data_recursive(
                        omsi_dependency_format=omsi_dependency_format,
                        omsi_main_parent=omsi_main_parent,
                        dependency_list=curr_full_dependency_list)
                    output_list = output_list + it_depend

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
        Get all direct and indirect dependencies associated with the analysis in form of a graph describing
        all nodes and links in the provenance hierarchy.

        :param include_omsi_dependency: Should the dependency_dict object be included in the entries
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
        :param prev_links: Previously established links in the list of nodes. Note, this list will be
                    modified by the call.
        :param parent_index: Index of the parent node in the list of prev_nodes for this call. May be
            None in case the parent we are calling this function for is not yet in the list. If None,
            then we will add our own parent that contains the dependencies to the list.
        :param metadata_generator: Optional parameter. Pass in a function that generates additional metadata
                    about a given omsi API object. Note, the key's level and path and name (i.e., name_key) are
                    already set by this function. The metadata_generator may overwrite these key's, however,
                    the path has to be unique as it is used to identify duplicate nodes. Overwriting the path
                    with a non-unique value, hence, will lead to errors (missing entries) when generating the graph.
                    Note, the metadata_generator function needs to support the following keyword arguments:

                          * in_dict : The dictionary to which the metadata should be added to.
                          * obj : The omsi file API object for which metadata should be generated
                          * name : A qualifying name for the object
                          * name_key : The key to be used for storing the name
        :param metadata_generator_kwargs: Dictionary of additional keyword arguments that should be passed to
                    the metadata_generator function.

        :returns: Dictionary containing two lists. 1) nodes : List of dictionaries, describing the elements
                  in the dependency graph. 2) links : List of tuples with the links in the graph. Each
                  tuple consists of two integer indices for the nodes list. For each node the following
                  entries are given:

                * dependency_dict: Optional key used to store the corresponding dependency_dict object.
                        Used only of include_omsi_dependency is True.
                * omsi_file_dependencydata: Optional key used to store the corresponding
                        omsi_file_dependencydata object. Used only of include_omsi_file_dependencydata is True.
                * name : Name of the dependency. The actual key is sepecified by name_key
                * level : The recursion level at which the object occurs.
                * path : The full path to the object
                * filename : The full path to the file
                * ... : Any other key/value pairs from the dependency_dict dict.

        """
        import os

        # 1) Initialize the nodes and the links
        if prev_nodes:
            nodes = prev_nodes
        else:
            nodes = []
        if prev_links:
            links = prev_links
        else:
            links = []

        # 2) Define internal helper functions
        # 2.1) Define a helper function to search the graph for a particular node
        def find_node_index(path_string, filename_string):
            """
            Internal helper function used to check whether a node already exists in the graph

            :param path_string: The path to object used for comparison.
            """
            # NOTE: The variable nodes is defined in the closure of the find_node_index(..) function, i.e.,
            #       nodes is defined at the beginning of the get_all_dependency_data_graph(...) that this
            #       function is contained in.
            if nodes:
                for i in range(len(nodes)):
                    if nodes[i]['path'] == path_string:
                        if omsi_file_common.same_file(nodes[i]['filename'], filename_string):
                            return i
            return None

        # 3) Add the original parent for which we are computing the graph if necessary
        if parent_index is None:
            omsi_main_parent = omsi_file_common.get_omsi_object(self.managed_group.parent)
            if not find_node_index(omsi_main_parent.name, omsi_main_parent.file.filename):
                nodes.append(self.__create_dependency_graph_node(
                    level=level,
                    name=os.path.basename(omsi_main_parent.name),
                    path=omsi_main_parent.name,
                    dependency_object=None,
                    omsi_object=omsi_main_parent,
                    include_omsi_dependency=include_omsi_dependency,
                    include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                    name_key=name_key,
                    metadata_generator=metadata_generator,
                    metadata_generator_kwargs=metadata_generator_kwargs))
                level += 1
                curr_index = len(nodes)-1
                parent_index = len(nodes)-1

        # 3) Iterate through all dependencies
        for item_obj in self.managed_group.items():
            try:
                # 3.1) Check if a node exists that represents the current dependency
                dep_obj = omsi_file_dependencydata(self.managed_group[unicode(item_obj[0])])
                omsi_obj = dep_obj.get_dependency_omsiobject(recursive=True)
                curr_path = omsi_obj.name
                curr_name = os.path.basename(curr_path)
                curr_index = find_node_index(omsi_obj.name, omsi_obj.file.filename)
                new_node = curr_index is None

                # 3.2) If no node exists for the current dependency then add a new node
                if new_node:
                    # 3.2.1)Add a new node with the basic metadata
                    nodes.append(self.__create_dependency_graph_node(
                        level=level,
                        name=curr_name,
                        path=curr_path,
                        dependency_object=dep_obj,
                        omsi_object=omsi_obj,
                        include_omsi_dependency=include_omsi_dependency,
                        include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                        name_key=name_key,
                        metadata_generator=metadata_generator,
                        metadata_generator_kwargs=metadata_generator_kwargs))
                    curr_index = len(nodes) - 1

                # 3) Add a link from the current node to its parent
                if parent_index is not None:
                    templink = {'source': parent_index,
                                'target': curr_index,
                                'link_type': dep_obj.get_dependency_type()}
                    # print templink
                    links.append(templink)

                # 4) If we can have recursive dependencies that we should trace,
                #    then trace them back and add them to the nodes and links lists
                #    NOTE: By checking for new_node we guarantee that the recursion
                #          is only performed the first time we encounter a particular node
                #    NOTE: By checking for the recursive parameter we ensure that recursion
                #          is only done if the user requests it
                if recursive and new_node:
                    # 4.1) We have a manager of dependencies and the manager has dependencies
                    if isinstance(omsi_obj, omsi_dependencies_manager) and omsi_obj.dependencies is not None:
                        # Recursively evaluate the dependencies to create the graph
                        nodes, links = omsi_obj.get_all_dependency_data_graph(
                            include_omsi_dependency=include_omsi_dependency,
                            include_omsi_file_dependencydata=include_omsi_file_dependencydata,
                            level=level + 1,
                            name_key=name_key,
                            prev_nodes=nodes,
                            prev_links=links,
                            parent_index=curr_index,
                            metadata_generator=metadata_generator,
                            metadata_generator_kwargs=metadata_generator_kwargs
                        )
            except:
                import sys
                warnings.warn("WARNING: Error occurred in omsi_file_dependencies::get_all_dependency_data_graph(...):" \
                      + str(sys.exc_info()))

        return nodes, links

    # 2.2) Define a helper function to create a new node for the graph
    @staticmethod
    def __create_dependency_graph_node(level,
                                       name,
                                       path,
                                       dependency_object,
                                       omsi_object,
                                       include_omsi_dependency=False,
                                       include_omsi_file_dependencydata=False,
                                       name_key='name',
                                       metadata_generator=None,
                                       metadata_generator_kwargs=None):
        """
        Internal helper function used to create a new node in the graph

        :param level: The recursion level at which the node exists
        :param name: The name of the node
        :param path: The path of the node
        :param dependency_object: The omsi_file_dependencydata object. May be None in case a node
            to a specific object is set
        :param omsi_object: The OpenMSI file API object. This is required and may NOT be None.
        :param include_omsi_dependency: Should the dependency_dict object be included in the entries
                    in the nodes dict?
        :param include_omsi_file_dependencydata: Should the omsi_file_dependencydata object be included in
                    the entries in the nodes dict?
        :param name_key: Which key should be used in the dicts to indicate the name of the object?
                    Default value is 'name'
        :param metadata_generator: Optional parameter. Pass in a function that generates additional metadata
                about a given omsi API object. Note, the key's level and path and name (i.e., name_key) are
                already set by this function. The metadata_generator may overwrite these key's, however,
                the path has to be unique as it is used to identify duplicate nodes. Overwriting the path
                with a non-unique value, hence, will lead to errors (missing entries) when generating the graph.
                Note, the metadata_generator function needs to support the following keyword arguments:

                      * in_dict : The dictionary to which the metadata should be added to.
                      * obj : The omsi file API object for which metadata should be generated
                      * name : A qualifying name for the object
                      * name_key : The key to be used for storing the name
        :param metadata_generator_kwargs: Dictionary of additional keyword arguments that should be passed to
                the metadata_generator function.

        :return: Dict describing the new node, containing the 'name', 'level', and 'path' and
            optionally 'dependency_dict' and/or 'omsi_file_dependencydata' and any additional
            data generated by the metadata_generator function
        """
        if metadata_generator_kwargs:
            metagenkwargs = metadata_generator_kwargs
        else:
            metagenkwargs = {}
        new_node = {'level': level,
                    name_key: name,
                    'path': path,
                    'filename': os.path.abspath(omsi_object.file.filename)}
        if include_omsi_dependency:
            if dependency_object:
                new_node['dependency_dict'] = dependency_object.get_omsi_dependency()
            else:
                new_node['dependency_dict'] = None
        if include_omsi_file_dependencydata:
            new_node['omsi_file_dependencydata'] = dependency_object
        #  Expand the metadata dict for the current object with  any user-defined metadata
        #  using the metadata generator function defined by the user
        if metadata_generator:
            metadata_generator(in_dict=new_node,
                               obj=omsi_object,
                               name=os.path.basename(omsi_object.name),
                               name_key=name_key,
                               **metagenkwargs)
        return new_node


##########################################################################
##########################################################################
#   The main class for accessing dependency datasets                     #
##########################################################################
##########################################################################
class omsi_file_dependencydata(omsi_file_common):

    """Class for managing data groups used for storing data dependencies"""
    @classmethod
    def __create__(cls, parent_group, dependency_data, use_relative_links=True):
        """
        Create a new dependency group within the given parent group and populate it
        with the required datasets using omsi_file_dependencydata.__populate_dependency__
        based on the given dependency information.

        :param parent_group: The h5py group opbject to which the dependency dataset should be added
        :type parent_group: h5py.Group
        :param dependency_data: The analysis dependency specification.
        :type dependency_data: omsi.shared.omsi_dependency_data
        :param use_relative_links: Should we use relative links to external files if possible. Default value
            is True.
        :type use_relative_links: bool

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
                                                                             dependency_data=dependency_data,
                                                                             use_relative_links=use_relative_links)
        return dependency_object

    @classmethod
    def __populate_dependency__(cls, dependency_group, dependency_data, use_relative_links=True):
        """
        Populate the given dependency_group with the associated data.

        :param dependency_group: The h5py group opbject to which the dependency dataset should be added
        :type dependency_group: h5py.Group
        :param dependency_data: The analysis dependency specification.
        :type dependency_data: omsi.shared.dependency_data
        :param use_relative_links: Should we use relative links to external files if possible. Default value
            is True.
        :type use_relative_links: bool

        :returns: omsi_file_dependencydata object for management of the newly created object

        """
        from omsi.datastructures.dependency_data import dependency_dict
        dep_group = dependency_group
        # 1) Save the name of the parameter
        param_name_data = dep_group.require_dataset(name=unicode(omsi_format_dependencydata.dependency_parameter),
                                                    shape=(1,),
                                                    dtype=omsi_format_common.str_type)
        help_attr = omsi_format_dependencydata.dependency_parameter_help_attr
        if omsi_format_common.str_type_unicode:
            param_name_data[0] = dependency_data['param_name']
            param_name_data.attrs[help_attr] = unicode(dependency_data['help'])
        else:
            param_name_data[0] = str(dependency_data['param_name'])
            param_name_data.attrs[help_attr] = str(dependency_data['help'])

        # 2) Save the selection
        selection_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_selection), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['selection'] is not None:
            from omsi.shared.data_selection import check_selection_string
            # This should always be True since dependency_dict checks for this
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

        # 3) Save the string describing the location of the main omsi object
        mainname_data = dep_group.require_dataset(name=unicode(omsi_format_dependencydata.dependency_mainname),
                                                  shape=(1,),
                                                  dtype=omsi_format_common.str_type)
        # 3.1) Get the name of the object
        mainname_string = dependency_data['omsi_object'].name
        # 3.2) Determine whether this is an internal dependency or an external file dependency
        mainname_filename = omsi_file_common.get_h5py_object(dependency_data['omsi_object']).file.filename
        if not omsi_file_common.same_file(mainname_filename, dep_group.file.filename):
            # 3.2.1 Create a relative path for the file if possible
            if use_relative_links:
                # Relative paths between the directories where the two files are stored
                relative_dir_path = os.path.relpath(os.path.dirname(os.path.abspath(mainname_filename)),
                                                    os.path.dirname(os.path.abspath(dep_group.file.filename)))
                # Create the relative path to the target file
                mainname_filename = os.path.join(relative_dir_path,
                                                 os.path.basename(mainname_filename))

            # 3.2.2 Create the string for the external link
            mainname_string = omsi_file_common.create_path_string(mainname_filename, mainname_string)

        # 3.3) Write the data to file
        mainname_string = unicode(mainname_string) if omsi_format_common.str_type_unicode else str(mainname_string)
        mainname_data[0] = unicode(mainname_string)

        # 4) Save the additional dataset name
        dataset_data = dep_group.require_dataset(name=unicode(
            omsi_format_dependencydata.dependency_datasetname), shape=(1,), dtype=omsi_format_common.str_type)
        if dependency_data['dataname']:
            if omsi_format_common.str_type_unicode:
                dataset_data[0] = unicode(dependency_data['dataname'])
            else:
                dataset_data[0] = str(dependency_data['dataname'])
        else:
            dataset_data[0] = u'' if omsi_format_common.str_type_unicode else ''

        # Save the dependency type if specified
        if dependency_data['dependency_type'] != dependency_dict.dependency_types['undefined']:
            dependency_type_data = dep_group.require_dataset(
                name=unicode(omsi_format_dependencydata.dependency_typename),
                shape=(1,),
                dtype=omsi_format_common.str_type)
            dependency_type_data[0] = dependency_data['dependency_type'] if \
                omsi_format_common.str_type_unicode else \
                str(dependency_data['dependency_type'])

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

    def get_dependency_objecttype(self, recursive=True):
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

    def get_parameter_help(self):
        """
        Get the help string for the parameter name if available.
        """
        try:
            group_name = unicode(omsi_format_dependencydata.dependency_parameter)
            attr_name = unicode(omsi_format_dependencydata.dependency_parameter_help_attr)
            return self.managed_group.attrs[attr_name]
        except:
            return ""

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

    def get_mainname(self):
        """
        Get the main name string describing the name of the object (and possibly path of the file if external)

        :return: String indicating the main name of the object that we link to
        """
        return self.managed_group[unicode(omsi_format_dependencydata.dependency_mainname)][0]

    def get_dependency_type(self):
        """
        Get the string describing the type of the dependency

        NOTE: If the type is missing in the file but we have a parameter name specified,
        then the default type 'parameter' will be returned other None is returned.

        :return: String indicating the type of the dependency or None if the type is not known.
        """
        if omsi_format_dependencydata.dependency_typename in self.managed_group.keys():
            return self.managed_group[unicode(omsi_format_dependencydata.dependency_typename)][0]
        else:
            from omsi.datastructures.dependency_data import dependency_dict
            param_name = self.get_parameter_name()
            if param_name is not None and len(param_name) > 0:
                return dependency_dict.dependency_types['parameter']
            return None

    def get_dependency_omsiobject(self, recursive=True, external_mode=None):
        """
        Get the omsi file API object corresponding to the object the dependency is pointing to.

        :param recursive: Should dependencies be resolved recursively, i.e., if the dependency points
                          to another dependencies. Default=True.
        :param external_mode: The file open mode (e.g., 'r', 'a') to be used when we encounter external dependencies,
            i.e., dependencies that are stored in external files. By default this is set to None,
            indicating that the same mode should be used in which this (i.e,. the current file
            describing the dependency) was opened. Allowed modes are 'r', 'r+', and 'a'. The modes
            'w', 'w+, 'x' are prohibited to ensure that we do not break external files.

        :returns: An omsi file API object (e.g., omsi_file_analysis or omsi_file_msidata) if the link
                  points to a group or the h5py.Dataset the link is pointing to.
        """
        # Get the name of the object
        dependency_path = self.get_mainname()
        # Determine if this is an internal or external dependency
        filename, omsi_object_name = omsi_file_common.parse_path_string(dependency_path)
        # Determine if we need to open a new file
        h5py_file = self.managed_group.file
        if filename is not None:
            if not omsi_file_common.same_file(filename, self.managed_group.file.filename):
                # Determine the filemode
                filemode = external_mode
                if filemode is None:
                    filemode = self.managed_group.file.mode
                filemode = filemode if filemode not in ['w', 'w+', 'x'] else 'a'
                # Resolve relative path names
                if not os.path.isabs(filename):
                    filename = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(self.managed_group.file.filename)),
                                                            filename))
                try:
                    h5py_file = h5py.File(filename, filemode)
                except IOError:
                    h5py_file = h5py.File(filename, 'r')
                    warnings.warn('External dependency forced to be opened in read-only mode.')

        h5py_object = h5py_file.file[unicode(omsi_object_name)]
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
        Get the dependency information as an omsi.shared.dependency_dict object
        (as defined in the omsi.shared.dependency_dict module)

        :returns: dependency_dict object with all the dependency data.
        """
        from omsi.datastructures.dependency_data import dependency_dict
        output_dependency = dependency_dict()
        try:
            output_dependency['param_name'] = self.get_parameter_name()
        except:
            warnings.warn('Failed to retrieve param_name for: ' + self.managed_group.name)
        try:
            output_dependency['link_name'] = self.get_link_name()
        except:
            warnings.warn('Failed to retrieve link_name for: ' + self.managed_group.name)
        try:
            output_dependency['selection'] = self.get_selection_string()
        except:
            warnings.warn('Failed to retrieve selection for: ' + self.managed_group.name)
        try:
            output_dependency['dependency_type'] = self.get_dependency_type()
        except:
            warnings.warn('Failed to retrieve dependency_type for: ' + self.managed_group.name)
        try:
            output_dependency['omsi_object'] = self.get_dependency_omsiobject()
        except:
            warnings.warn('Failed to retrieve omsi_object for: ' + self.managed_group.name)
        try:
            output_dependency['dataname'] = self.get_dataset_name()
        except:
            warnings.warn('Failed to retrieve dataname for: ' + self.managed_group.name)
        try:
            output_dependency['help'] = self.get_parameter_help()
        except:
            warnings.warn('Failed to retrieve help for: ' + self.managed_group.name)
        return output_dependency
