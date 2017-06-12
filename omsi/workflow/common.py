"""
Module defining basic data structures used by workflows.
"""
from omsi.shared.log import log_helper
import argparse


class RawDescriptionDefaultHelpArgParseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                                 argparse.RawDescriptionHelpFormatter):
    """
    Simple derived formatter class for use with argparse used by the
    cl_analysis_driver class. This formatter combines the default
    argparse.ArgumentDefaultsHelpFormatter and
    argparse.RawDescriptionHelpFormatter
    for formatting arguments and help descriptions.

    """
    pass


class analysis_task_list(list):
    """
    Define a python list of analyses to be executed as a workflow. The list allows only for
    storage of analysis_base objects and ensures uniqueness of elements in the list.

    """
    @classmethod
    def from_script_files(cls,
                          script_files):
        """
        Same as from_script, but the script is read from the given files.

        :param script_files: List of strings with the paths to the script files. If only a single
            script is used, then a single string may be used as well.

        :return: An instance of workflow_base with the specification of the workflow to be executed
        """
        scripts = []
        inscripts = [script_files, ] if isinstance(script_files, basestring) else script_files
        for script_path in inscripts:
            infile = open(script_path, 'r')
            scripttext = infile.read()
            scripts.append(scripttext)
            infile.close()
        return cls.from_scripts(scripts=scripts)

    @classmethod
    def from_scripts(cls,
                     scripts):
        """
        Evaluate the workflow script to extract all analyses to be run.

        NOTE: This function executes using eval(..), i.e., there are NO safeguards against malicious codes.

        :param scripts: The script with the setup of the workflow. This should only include
            the definition of analyses and their inputs.

        :return: An instance of workflow_base with the specification of the workflow to be executed
        """

        # Check which analyses exist already as we only want to include the analyses created by the script
        from omsi.analysis.base import analysis_base
        current_analyses = cls.all()
        # Evaluate the workflows script
        if isinstance(scripts, basestring):
            exec(scripts)
        else:
            for script in scripts:
                exec(script)
        # Create the list of new analyses that have been created
        new_analysis_task_list = cls()
        for ana_obj in analysis_base.get_analysis_instances():
            if ana_obj not in current_analyses:
                new_analysis_task_list.add(ana_obj)
        return new_analysis_task_list

    @classmethod
    def all(cls):
        """
        Get an analysis_task_list of all current analysis_base objects

        :return: New analysis_task_list of all current analysis objects
        """
        from omsi.analysis.base import analysis_base
        current_analyses = analysis_task_list()
        for ana_obj in analysis_base.get_analysis_instances():
            current_analyses.append(ana_obj)
        return current_analyses

    def __init__(self,
                 analysis_objects=None):
        """
        Initialize the set of analysis tasks to be performed as part of the workflow

        :param analysis_objects: Set or list of unique analysis objects to be added.
            Duplicates will be removed.
        """
        super(analysis_task_list, self).__init__()
        self.update(analysis_objects)

    def __getitem__(self, item):
        """
        Get the analysis object with the given integer index of analysis_identifier
        :param item: Integer index of the analysis or string with the analysis identifier to be retrieved.
            NOTE: If analyses identifiers are user-defined, and if they are not unique, then the subset
            with matching identifiers will be retrieved.

        :return: Instance of omsi.analysis.base.analysis_base with the matching analysis object. In case that
            multiple analyses match the key, then a analysis_task_list object will be returned instead
        """
        if isinstance(item, int):
            return list(self)[item]
        elif isinstance(item, basestring):
            matching_analyses = []
            for ana_obj in self:
                if ana_obj.analysis_identifier == item:
                    matching_analyses.append(ana_obj)
            if len(matching_analyses) == 1:
                return matching_analyses[0]
            if len(matching_analyses) > 1:
                return analysis_task_list(matching_analyses)
        else:
            raise KeyError('Invalid key type')
        raise KeyError('Invalid key')

    def update(self, analysis_objects):
        """
        Return the set with elements added from the given set of analysis_objects.

        This is the same as set.update() but we ensure that only analysis_base objects
        are added.

        :param analysis_objects: List or set of analysis_base objects to be added to workflow set

        :raise: ValueError is raised in case that objects that are not instances of analysis_base are to be added.

        :return: self with elements added to self.
        """
        from omsi.analysis.base import analysis_base
        # Nothing to be added
        if analysis_objects is None:
            return self
        # Check that all elements to be added are instances of analysis_base
        for ana_obj in analysis_objects:
            self.append(ana_obj)
        return self

    def add(self, analysis_object):
        """
        Same as append
        """
        self.append(analysis_object)

    def clear(self):
        """
        Remove all elements from the list
        """
        self[:] = []

    def append(self, analysis_object):
        """
        Add a given analysis to the set of object to be executed by the workflow

        This is the same as set.add() but we ensure that only analysis_base objects
        are added.

        :param analysis_object: Analysis object to be added to the execution.
            All dependencies of the analysis will also be executed as part of the
            execution.
        :type analysis_object: omsi.analysis.base.analysis_base

        :raises: ValueError is raised if the given analysis_object is invalid
        """
        from omsi.analysis.base import analysis_base
        if isinstance(analysis_object, analysis_base):
            if analysis_object in self:
                log_helper.debug(__name__, "Analysis already in the list of tasks")
                return
            log_helper.debug(__name__, "Adding analysis object to the workflow set. " + str(analysis_object))
            super(analysis_task_list, self).append(analysis_object)
        else:
            raise ValueError('Analysis is not of type omsi.analysis.base.analysis_base')

    def insert(self, index, analysis_object):
        """
        Insert a given analysis object at the given location

        :param index: Location where the obejct should be inserted
        :param analysis_object: The analysis object to be inserted

        """
        from omsi.analysis.base import analysis_base
        if isinstance(analysis_object, analysis_base):
            if analysis_object in self:
                log_helper.debug(__name__, "Analysis already in the list of tasks")
                return
            log_helper.info(__name__, "Inserting analysis object in the workflow list. " + str(analysis_object))
            super(analysis_task_list, self).insert(index, analysis_object)
        else:
            raise ValueError('Analysis is not of type omsi.analysis.base.analysis_base')

    def add_analysis_dependencies(self):
        """
        Add the dependencies of all analyses to the workflow list in case they are missing.

        This function is recursive, step-by-step adding all dependencies of the workflow to the list of
        tasks to be executed, until no more dependencies are found.

        Usually this function is called by the workflow executor itself before running the analysis and
        should not need to be called by the user.

        :returns: Integer indicating the number of dependencies added to the list of tasks
        """
        # Add all direct dependencies
        add_objects = self.get_additional_analysis_dependencies()
        for new_analysis in add_objects:
            self.add(new_analysis)

        # Recursively add all indirect dependencies if add new objects
        added_dependencies = len(add_objects)
        if added_dependencies > 0:
            added_dependencies += self.add_analysis_dependencies()
        # Return the number of dependencies we added
        return added_dependencies

    def get_additional_analysis_dependencies(self):
        """
        Compute a list of all dependencies of the current list of analyses (excluding analyses that
        are already in the the list of tasks.

        :return: analysis_task_list of all analysis dependencies
        """
        from omsi.dataformat.omsi_file.common import omsi_file_common
        from omsi.analysis.base import analysis_base

        missing_dependencies = analysis_task_list()
        for analysis_obj in self:
            for dependency_param_obj in analysis_obj.get_all_dependency_data():
                dependency_analysis = dependency_param_obj['data']['omsi_object']
                if isinstance(dependency_analysis, analysis_base):
                    if dependency_analysis not in self:
                        missing_dependencies.add(dependency_analysis)
                elif isinstance(dependency_analysis, omsi_file_common):
                    pass  # Ignore dependencies on data files. We do not need to execute those
                else:
                    log_helper.warning(__name__, 'Unknown dependency object type that cannot be processed by workflow.'
                                       + str(dependency_param_obj))
        return missing_dependencies

    def get_all_dependency_data(self):
        """
        Get the complete list of all direct and indirect dependencies of all analysis tasks.

        NOTE: These are only the direct dependencies as specified by the analysis itself.
        Use  get_all_dependency_data_recursive(..) to also get the indirect dependencies of
        the analysis due to dependencies of the dependencies themselves.

        :returns: List of parameter_data objects that define dependencies.

        """
        all_analyses = self.get_additional_analysis_dependencies()
        all_analyses.update(self)
        dependencies = [ana.get_all_dependency_data() for ana in all_analyses]
        return dependencies

    def get_all_parameter_data(self,
                               exclude_dependencies=False):
        """
        Get the complete list of all parameters

        :param exclude_dependencies: Boolean indicating whether we should exclude parameters
            that define dependencies from the list

        :return: List of omsi.analysis.analysis_data.parameter_data objects with the description of the parameters
        """
        parameters = [ana.get_all_parameters(exclude_dependencies) for ana in self]
        return parameters

    def get_all_run_info(self):
        """
        Get a list of dict with the complete info about the last run of each of the analysis analysis

        :return: List of run_info_dict objects, one for each analysis

        """
        run_infos = [ana.get_all_run_info() for ana in self]
        return run_infos

    def get_all_analysis_data(self):
        """
        Get a list of all output data objects for all analysis

        :return: List of omsi.analysis.analysis_data.analysis_data objects, one for each analysis
        """
        analysis_data = [ana.get_all_analysis_data() for ana in self]
        return analysis_data

    def get_all_analysis_identifiers(self):
        """
        Get a list of all analysis identifiers

        :return: List of strings with the analysis identifier of each analysis
        """
        identifiers = [ana.get_analysis_identifier() for ana in self]
        return identifiers

    def add_all(self):
        """
        Add all known analyses to the workflow list

        :return: The updated analysis_task_list with all analyses added
        """
        from omsi.analysis.base import analysis_base
        for ana in analysis_base.get_analysis_instances():
            self.add(ana)
        return self

    def enable_time_and_usage_profiling(self, enable=True):
        """
        Enable or disable profiling of time and usage of code parts of execute_analysis for all analyses.

        :param enable: Enable (True) or disable (False) profiling
        :type enable: bool

        :raises: ImportError is raised if a required package for profiling is not available.
        """
        for ana in self:
            ana.enable_time_and_usage_profiling(enable=enable)

    def enable_memory_profiling(self, enable=True):
        """
        Enable or disable line-by-line profiling of memory usage of execute_analysis.

        :param enable_memory: Enable (True) or disable (False) line-by-line profiling of memory usage
        :type enable_memory: bool

        :raises: ImportError is raised if a required package for profiling is not available.
        """
        for ana in self:
            ana.enable_memory_profiling(enable=enable)

    def set_undefined_analysis_identifiers(self):
        """
        Check that all analysis descriptors are set to a value different than "undefined" and
        set the descriptor based on their index in the list if necessary.
        """
        ana_index = 0
        for ana in self:
            if not ana.analysis_identifier_defined():
                ana.set_analysis_identifier('ana_%i' % ana_index)
            ana_index += 1

    def analysis_identifiers_unique(self):
        """
        Check whether all identifiers of the analyses in the this list are unique.
        :return: bool
        """
        identifiers = self.get_all_analysis_identifiers()
        unique_identifiers = list(set(identifiers))
        identifiers_unique = len(identifiers) == len(unique_identifiers)
        return identifiers_unique

    def make_analysis_identifiers_unique(self):
        """
        Update analysis identifiers to be unique.

        Side effects: This function updates the analysis tasks stored in the set

        :return: self, i.e., the modified object with identifiers updated
        """
        identifiers = self.get_all_analysis_identifiers()
        unique_identifiers = list(set(identifiers))
        num_update = len(identifiers) - len(unique_identifiers)
        if num_update > 0:
            log_helper.debug(__name__, "%i analyses have non-unique identifiers and will be updated" % num_update)
        ana_index = 0
        for ana in self:
            current_identifier = ana.get_analysis_identifier()
            if current_identifier not in unique_identifiers:
                ana.set_analysis_identifier('ana_' + str(ana_index) + "_" + unicode(current_identifier))
            ana_index += 1
        return self

    def task_status_stats(self):
        """
        Compute the number of tasks that are complete, ready to run, and waiting on dependencies
        :return: Tuple of ints with:
            * number of completed tasks
            * total number of tasks waiting to be executed
            * total number of waiting tasks that are ready to be run
            * total number of waiting tasks that are blocked, e.g., cannot run due to unresolved dependencies

        """
        waiting_tasks = 0
        ready_to_run_tasks = 0
        blocked_tasks = 0
        completed_tasks = 0
        for analysis in self:
            if analysis.update_analysis:
                waiting_tasks += 1
                if len(analysis.check_ready_to_execute()) == 0:
                    ready_to_run_tasks += 1
                else:
                    blocked_tasks += 1
            else:
                completed_tasks += 1
        return completed_tasks, waiting_tasks, ready_to_run_tasks, blocked_tasks

    def get_blocking_tasks(self):
        """
        Get the tasks that are blocking other from running
        :return: List of blocking tasks
        """
        from omsi.analysis.base import analysis_base
        blocking_tasks = []
        for analysis in self:
            if analysis.update_analysis:
                pending_inputs = analysis.check_ready_to_execute()
                for pi in pending_inputs:
                    if isinstance(pi['data']['omsi_object'], analysis_base):
                        blocking_tasks.append(pi['data']['omsi_object'])
        blocking_tasks = list(set(blocking_tasks))
        return blocking_tasks

