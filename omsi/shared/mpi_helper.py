try:
    from mpi4py import MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False
import numpy as np
import itertools
import warnings
import time
import os


class parallel_over_axes(object):
    """
    Helper class used to parallelize the execution of a function using MPI by splitting the
    input data into sub-blocks along a given set of axes.

    :ivar task_function: The function we should run.
    :ivar task_function_params: Dict with the input parameters for the function.
        may be None or {} if no parameters are needed.
    :ivar main_data: Dataset over which we should parallelize
    :ivar split_axes:  List of integer axis indicies over which we should parallelize
    :ivar main_data_param_name: The name of data input parameter of the task function
    :ivar root: The master MPI rank (Default=0)
    :ivar schedule: The task scheduling schema to be used (see parallel_over_axes.SCHEDULES
    :ivar collect_output: Should we collect all the output from the ranks on the master rank?
    :ivar schedule: The parallelization schedule to be used. See also parallel_over_axes.schedule
    :ivar result: The result form the task_function. If collect_output is set and we are the root
        then this will be the output of all tasks
    :ivar comm: The MPI communicator used for the parallelization. Default value is MPI.COMM_WORLD

    """
    SCHEDULES = {'STATIC_1D': 'STATIC_1D',
                 'STATIC': 'STATIC',
                 'DYNAMIC': 'DYNAMIC'}

    MPI_MESSAGE_TAGS = {'RANK_MSG': 11,
                        'BLOCK_MSG': 12,
                        'COLLECT_MSG': 13}

    def __init__(self,
                 task_function,
                 task_function_params,
                 main_data,
                 split_axes,
                 main_data_param_name,
                 schedule=SCHEDULES['STATIC_1D'],
                 root=0,
                 comm=None):
        """

        :param task_function: The function we should run.
        :param task_function_params: Dict with the input parameters for the function.
            May be None or {} if no parameters are needed.
        :param main_data: Dataset over which we should parallelize
        :param split_axes:  List of integer axis indicies over which we should parallelize
        :param main_data_param_name: The name of data input parameter of the task function
        :param root: The master MPI rank (Default=0)
        :param schedule: The task scheduling schema to be used (see parallel_over_axes.SCHEDULES
        :param comm: The MPI communicator used for the parallelization. Default value is None, in which case
            MPI.COMM_WORLD is used

        """
        if not is_mpi_available():
            raise ValueError("MPI is not available. MPI is required for parallel execution.")
        self.task_function = task_function
        self.schedule = schedule
        self.split_axes = split_axes
        self.main_data = main_data
        self.main_data_param_name = main_data_param_name
        self.task_function_params = task_function_params
        if self.task_function_params is None:
            self.task_function_params = {}
        self.root = root
        self.result = None
        self.blocks = None
        self.comm = get_comm_world() if comm is None else comm

    def run(self):
        """
        Call this function to run the function in parallel.

        :return: Tuple with the following elements:

            1) List of tuples with the results from the local execution of the task_function. Each
              tuple is the result from one return of the task_function. In the case of static
              execution, this is always a list of length 1.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition whereas in the
               case of dynamic scheduling we usually have single integer point selections
               for each task.

        """
        if self.schedule == self.SCHEDULES['DYNAMIC']:
            return self.__run_dynamic()
        elif self.schedule == self.SCHEDULES['STATIC_1D'] or self.schedule == self.SCHEDULES['STATIC']:
            return self.__run_static_1D()

    def collect_data(self):
        """
        Collect the results from the parallel execution to the self.root rank.

        :return: On worker ranks (i.e., MPI_RANK!=self.root) this is simply the
            self.result and self.blocks containing the result created by run function.
            On the root rank (i.e., MPI_RANK!=self.root) this is a tuple of two lists
            containing with the self.result and self.blocks from all ranks respectively.

        """
         # Collect the output
        rank = get_rank()
        start_time = time.time()
        if rank == self.root:
            print "COLLECTING RESULTS"
        collected_data = self.comm.gather(self.result, root=self.root)
        collected_blocks = self.comm.gather(self.blocks, root=self.root)
        if rank == self.root:
            self.result = collected_data
            self.blocks = collected_blocks
        if self.schedule == self.SCHEDULES['DYNAMIC']:
            temp = list(self.result)
            temp.pop(self.root)
            self.result = temp
            temp = list(self.blocks)
            temp.pop(self.root)
            self.blocks = temp
        end_time = time.time()
        run_time = end_time - start_time
        if rank == self.root:
            print "TIME FOR COLLECTING DATA FROM ALL TASKS: " + str(run_time)
        return self.result, self.blocks

    def __run_static_1D(self):
        """
        Run the task function using a static task decomposition schema.

        The data is divided into sub-blocks along the largest split_axis

        :return: Tuple with the following elements:

            1) List of tuples with the results from the local execution of the task_function. Each
              tuple is the result from one return of the task_function. In the case of static
              execution, this is always a list of lenght 1.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition.

        """
        start_time = time.time()
        # Get MPI parameters
        rank = get_rank()
        size = get_size()

        # Get data shape parameters and compute the data blocks
        # Determine the longest axis along which we can split the data
        axes_shapes = np.asarray(self.main_data.shape)[self.split_axes]
        total_num_subblocks = np.prod(axes_shapes)
        if total_num_subblocks < size:
            size = total_num_subblocks
            if rank == self.root:
                print "Insufficient number of blocks for number of MPI ranks. Some ranks will remain idle"
        axes_sort_index = np.argsort(axes_shapes)[::-1]
        split_axis = self.split_axes[axes_sort_index[0]]
        split_axis_size = axes_shapes[split_axis]
        if split_axis_size < size:
            raise NotImplementedError("STATIC scheduling currently parallelizes only over one axis, " +
                                      "and the largest axis is too small to fill all MPI tasks")
        # Determine the size of 1D block
        block_size = int(split_axis_size / float(size) + 0.5)
        if block_size * size > split_axis_size and block_size > 1:
            block_size -= 1

        # Compute a block for every rank
        self.blocks = [slice(None)] * len(self.main_data.shape)
        start_index = rank * block_size
        stop_index = start_index + block_size
        if rank == (size-1):
            if stop_index != split_axis_size:
                stop_index = split_axis_size
        self.blocks[axes_sort_index[0]] = slice(start_index, stop_index)
        self.blocks = tuple(self.blocks)
        print "Rank: " + str(rank) + " Block: " + str(self.blocks)

        # Execute the task_function on the given data block
        task_params = self.task_function_params
        task_params[self.main_data_param_name] = self.main_data[self.blocks]
        self.result = self.task_function(**task_params)

        end_time = time.time()
        run_time = end_time - start_time
        print "TIME FOR PROCESSING THE DATA BLOCK: " + str(run_time)

        # Return the output
        self.result = [self.result,]
        self.blocks = [self.blocks,]
        return self.result, self.blocks

    def __run_dynamic(self):
        """
        Run the task function using dynamic task scheduling.

        The root rank divides the data into sub-tasks and sends the tasks to available MPI
        processes on request.

        :return: Tuple with the following elements:

            1) List of tuples with the results from the local execution of the task_function. Each
              tuple is the result from one return of the task_function.
            2) List of block_indexes. Each block_index is a tuple with the selection used to
               divide the data into sub-blocks. In the case of static decomposition we have
               a range slice object along the axes used for decomposition.

        """
        import time
        rank = get_rank()
        size = get_size()

        if size < 2:
            warnings.warn('DYNAMIC task scheduling requires at least 2 MPI ranks. Using STATIC scheduling instead.')
            return self.__run_static_1D()

        # We are the controlling rank
        if rank == self.root:
            self.result = []
            self.blocks = []
            # Get data shape parameters and compute the data blocks
            axes_shapes = np.asarray(self.main_data.shape)[self.split_axes]
            total_num_subblocks = np.prod(axes_shapes)
            if total_num_subblocks < size:
                size = total_num_subblocks
                if rank == self.root:
                    warnings.warn("Insufficient number of blocks for number of MPI ranks. Some ranks will remain idle")

            # Compute the list of all possible blocks
            base_blocks = [[slice(None)]] * len(self.main_data.shape)
            for axis_index in self.split_axes:
                base_blocks[axis_index] = range(self.main_data.shape[axis_index])
            block_tuples = itertools.product(*base_blocks)

            # Communicate blocks with task ranks
            print "PROCESSING DATA BLOCKS"
            start_time = time.time()
            block_index = 0
            for block_selection in block_tuples:
                request_rank = self.comm.recv(source=MPI.ANY_SOURCE, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                self.comm.send((block_index, block_selection),
                               dest=request_rank,
                               tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                block_index += 1
                if (block_index % 100) == 0:
                    print (block_index, total_num_subblocks, request_rank)
            end_time = time.time()
            run_time = end_time - start_time
            print "TIME FOR RUNNING TASK FUNCTION: " + str(run_time)
            start_time = time.time()
            print "FINALIZING"
            # Terminate all ranks and receive all data from the different ranks if requested
            all_ranks_status = np.zeros(size, 'bool')
            all_ranks_status[self.root] = True
            while not np.all(all_ranks_status):

                request_rank = self.comm.recv(source=MPI.ANY_SOURCE, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                self.comm.send((None, None), dest=request_rank, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                all_ranks_status[request_rank] = True

            end_time = time.time()
            run_time = end_time - start_time
            print "TIME FOR FINALIZING TASKS: " + str(run_time)

        # We are a rank that has to run tasks
        else:
            # Request a new data block
            self.result = []
            self.blocks = []

            while True:
                self.comm.send(rank, dest=self.root, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                block_index, block_selection = self.comm.recv(source=self.root, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                if block_index is None:
                    break
                # Execute the task_function on the given data block
                task_params = self.task_function_params
                task_params[self.main_data_param_name] = self.main_data[block_selection]
                self.result.append(self.task_function(**task_params))
                self.blocks.append(block_selection)

        # Return the result
        return self.result, self.blocks



def imports_mpi(python_object):
    """
    Check whether the given class import mpi

    The implementation inspects the source code of the
    analysis to see if MPI is imported by the code.
    """
    import inspect
    import re
    code = inspect.getsource(python_object)
    object_imports_mpi = re.search('import\s+mpi4py', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('from\s+mpi4py\s+import', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('import\s+mpi4py', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('import\s+mpi', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('from\s+mpi\s+import', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('from\s+omsi.shared.mpi_helper\s+import', code) is not  None
    object_imports_mpi = object_imports_mpi or re.search('import\s+omsi.shared.mpi_helper', code) is not  None
    return object_imports_mpi


def gather(data, comm=None, root=0):
    """
    MPI gather operation or return a list with just [data,] if MPI is not available

    :param data: The data to be gathered
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :param root: The rank where the data should be collected to. Default value is 0
    :return: List of data objects from all the ranks
    """
    if MPI_AVAILABLE:
        my_comm = comm if comm is not None else MPI.COMM_WORLD
        return comm.gather(data, root)
    else:
        return [data, ]

def barrier(comm=None):
    """
    MPI barrier operation or no-op when running without MPI

    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    """
    if MPI_AVAILABLE:
        my_comm = comm if comm is not None else MPI.COMM_WORLD
        my_comm.barrier()
    else:
        pass

def get_rank(comm=None):
    """
    Get the current process rank
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :return: The integer index of the rank
    """
    rank = 0
    if MPI_AVAILABLE:
        if comm:
            rank = comm.Get_rank()
        else:
            rank = MPI.COMM_WORLD.Get_rank()
    return rank


def get_size(comm=None):
    """
    Get the size of the current communication domain/
    :param comm: MPI communicator. If None, then MPI.COMM_WORLD will be used.
    :return: The integer index of the rank
    """
    size = 1
    if MPI_AVAILABLE:
        if comm:
            size = comm.Get_size()
        else:
            size = MPI.COMM_WORLD.Get_size()
    return size


def get_comm_world():
    """
    Get MPI.COMM_WORLD
    :return: mpi communicator or None if MPI is not available
    """
    if MPI_AVAILABLE:
        return MPI.COMM_WORLD
    else:
        return None


def is_mpi_available():
    """
    Check if MPI is available. Same as MPI_AVAILABLE
    :return: bool indicating whether MPI is available
    """
    return MPI_AVAILABLE

def mpi_type_from_dtype(dtype):
    """
    Ge the the corresponding MPI type for the given basic numpy dtype

    :param dtype: Basic numpy dtype to be mapped to the MPI type
    :return: The MPI type or None if not found
    """
    if MPI_AVAILABLE:
        try:
            return MPI.__TypeDict__[dtype.char]
        except AttributeError:  # Older versions of mpi4py use _typedict
            return MPI._typedict[dtype.char]
        except KeyError:
            return None
        except:
            return None
    else:
        return None


# class suppress_stdout_and_stderr(object):
#     """
#     Simple context manager that can be used to suppress the output to stdout and stderr.
#
#     To suppress the output for a set of calls use:
#
#     >>> with suppress_stdout_and_stderr():
#     >>>     my_function(...)
#
#     """
#
#     def __init__(self):
#         # Open a pair of null files
#         self.null_out_files = {'stdout': os.open(os.devnull,os.O_RDWR),
#                                'stderr': os.open(os.devnull,os.O_RDWR)}
#         # Save the original output streams
#         self.orginal_out_files = {'stdout': os.dup(1),
#                                   'stderr': os.dup(2)}
#
#     def __enter__(self):
#         # Suppress the output to stdout and stderr by assigning null-fiules to the outputs
#         os.dup2(self.null_out_files['stdout'], 1)
#         os.dup2(self.null_out_files['stderr'], 2)
#
#     def __exit__(self, *_):
#         # Restore the stdout and stderr to the original function
#         os.dup2(self.original_out_files['stdout'], 1)
#         os.dup2(self.original_out_files['stderr'], 2)
#         # Close the null files
#         for nullfile in self.null_out_files.values():
#             os.close(nullfile)
#         self.null_out_files = {'stdout': None, 'stderr': None}



