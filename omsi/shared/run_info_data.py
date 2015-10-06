"""
Module with helper data structures for recording runtime provenance data
"""

import warnings
import platform
import datetime
import sys
import omsi.shared.mpi_helper as mpi_helper
import time
from omsi.shared.log import log_helper


class run_info_dict(dict):
    """
    Simple dictionary class for collecting runtime information
    """
    DEFAULT_TIME_FORMAT = '%Y-%m-%d %H:%M:%S.%f'

    def __init__(self, *args, **kwargs):
        super(run_info_dict, self).__init__(*args, **kwargs)

    def record_preexecute(self,
                          root=0,
                          comm=None):
        """
        Record basic runtime information in this dict before the exeuction is started.


        Function used to record runtime information prior to executing the process we want to track, e.g.,
        the `execute_analysis(...)` of a standard analysis.

        The function may be overwritten in child classes to add recording of
        additional runtime information. All runtime data should be recorded in the
        main dict (i.e, self). This ensures in the case of standard analysis that
        the data is stored in the HDF5 file. Other data should be stored in separate
        variables that we may add to the object.

        When overwriting the function we should typically call super(...,self).runinfo_record_pretexecute()
        last in the custom version to ensure that the start_time is properly recorded right before
        the execution of the analysis.

        :param root: Used for logging only. The rank where logging should be performed.
            Default is None in which case logging is performed on all ranks.

        :param comm: Used for logging only. The MPI communicator to be used. Default value is None,
            in which case MPI.COMM_WORLD is used.

        """
        log_helper.debug(__name__, 'Recording pre-execution runtime data', root=root, comm=comm)
        # Record basic runtime environment information using the platform module
        try:
            self['architecture'] = unicode(platform.architecture())
            self['java_ver'] = unicode(platform.java_ver())
            self['libc_ver'] = unicode(platform.libc_ver())
            self['linux_distribution'] = unicode(platform.linux_distribution())
            self['mac_ver'] = unicode(platform.mac_ver())
            self['machine'] = unicode(platform.machine())
            self['node'] = unicode(platform.node())
            self['platform'] = unicode(platform.platform())
            self['processor'] = unicode(platform.processor())
            self['python_branch'] = unicode(platform.python_branch())
            self['python_build'] = unicode(platform.python_build())
            self['python_compiler'] = unicode(platform.python_compiler())
            self['python_implementation'] = unicode(platform.python_implementation())
            self['python_revision'] = unicode(platform.python_revision())
            self['python_version'] = unicode(platform.python_version())
            self['release'] = unicode(platform.release())
            self['system'] = unicode(platform.system())
            self['uname'] = unicode(platform.uname())
            self['version'] = unicode(platform.version())
            self['win32_ver'] = unicode(platform.win32_ver())
        except:
            warnings.warn("WARNING: Recording of execution provenance failed: " + str(sys.exc_info()))

        # Attempt to record the svn version information
        try:
            import subprocess
            self['svn_ver'] = subprocess.check_output('svnversion').rstrip('\n')
        except:
            warnings.warn("Recording of svn version information failed: "+str(sys.exc_info()))

        # Record the start time for the analysis
        self['start_time'] = unicode(datetime.datetime.now())

    def record_postexecute(self,
                           execution_time=None,
                           root=0,
                           comm=None):
        """
        Function used to record runtime information after the task we want to track is comleted, e.g.
        the `execute_analysis(...)` function of a standard analysis.

        The function may be overwritten in child classes to add recording of
        additional runtime information.

        When overwriting the function we should call super(...,self).runinfo_record_postexecute(execution_time)
        in the custom version to ensure that the execution and end_time are properly
        recorded.

        :param execution_time: The total time it took to execute the analysis. May be None, in which
            case the function will attempt to compute the execution time based on the start_time
            (if available) and the the current time.

        :param root: Used for logging only. The rank where logging should be performed.
            Default is None in which case logging is performed on all ranks.

        :param comm: Used for logging only. The MPI communicator to be used. Default value is None,
            in which case MPI.COMM_WORLD is used.

        """
        log_helper.debug(__name__, 'Recording post-execution runtime data', root=root, comm=comm)
        # Finalize recording of post execution provenance
        self['end_time'] = unicode(datetime.datetime.now())
        if execution_time:
            self['execution_time'] = unicode(execution_time)
        elif 'start_time' in self:
            start_time = run_info_dict.string_to_time(self['start_time'])
            stop_time = run_info_dict.string_to_time(self['end_time'])
            execution_time = stop_time - start_time    # TODO: This only gives execution time in full seconds right now
        else:
            self['execution_time'] = None

    def clean_up(self,
                 root=0,
                 comm=None):
        """
        Clean up the runinfo object. In particular remove empty keys that
        either recorded None or recorded just an empty string.

        This function may be overwritten to also do clean-up needed
        due to additional custom runtime instrumentation.

        When overwriting this function we should call super(..., self).runinfo_clean_up()
        at the end of the function to ensure that the runinfo dictionary
        is clean, i.e., does not contain any empty entries.

        :param root: Used for logging only. The rank where logging should be performed.
            Default is None in which case logging is performed on all ranks.

        :param comm: Used for logging only. The MPI communicator to be used. Default value is None,
            in which case MPI.COMM_WORLD is used.
        """
        log_helper.debug(__name__, 'Clean up runtime data', root=root, comm=comm)
        # Remove empty items from the run_info dict
        for ri_key, ri_value in self.items():
            try:
                if ri_value is None or len(ri_value) == 0:
                    self.pop(ri_key)
            except:
                pass

    def gather(self,
               root=0,
               comm=None):
        """
        Simple helper function to gather the runtime information---that has been collected on
        multiple processes when running using MPI---on a single root process

        :param root: The process where the runtime information should be collected.
            Default is None in which case 0 is used

        :param comm: The MPI communicator to be used. Default value is None,
            in which case MPI.COMM_WORLD is used..

        :return: If we have more than one processes then this function returns a
            dictionary with the same keys as usual for the run_info but the
            values are now lists with one entry per mpi processes. If we only have
            a single process, then the run_info object will be returned without
            changes. NOTE: Similar to mpi gather, the function only collects
            information on the root. All other processes will return just their
            own private runtime information.

        """
        if mpi_helper.MPI_AVAILABLE:
            if not comm:
                comm = mpi_helper.get_comm_world()
            if comm.Get_size() > 1:
                log_helper.debug(__name__, 'Gather runtime data from parallel tasks', root=root, comm=comm)
                self['mpi_rank'] = comm.Get_rank()
                run_data = comm.gather(self, root=root)
                if comm.Get_rank() == root:
                    merged_run_data = {}
                    for run_dict in run_data:
                        for key in run_dict:
                            try:
                                merged_run_data[key].append(run_dict[key])
                            except KeyError:
                                merged_run_data[key] = [run_dict[key]]
                    return merged_run_data
        return self

    @staticmethod
    def string_to_structime(time_string, time_format=None):
        """
        Covert a time string to a time.struct_time using time.strptime

        :param time_string: String with the time, e.g, with the start time of a program.
        :param time_format: The time format to be used or None in which case run_info_dict.DEFAULT_TIME_FORMAT
            will be used.

        """
        return time.strptime(time_string,
                             time_format if time_format is not None else run_info_dict.DEFAULT_TIME_FORMAT)

    @staticmethod
    def string_to_time(time_string, time_format=None):
        """
        Convert a time string to local time object using time.mktime.

        :param time_string: String with the time, e.g, with the start time of a program.
        :param time_format: The time format to be used or None in which case run_info_dict.DEFAULT_TIME_FORMAT
            will be used.

        """
        return time.mktime(time.strptime(time_string,
                                         time_format if time_format is not None else run_info_dict.DEFAULT_TIME_FORMAT))

