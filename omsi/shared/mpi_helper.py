try:
    from mpi4py import MPI
    MPI_AVAILABLE = True
except ImportError:
    MPI_AVAILABLE = False
import numpy as np
import itertools
import warnings
import sys

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
                 'STATIC_2D': 'STATIC_2D',
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
                 collect_output=True,
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
        :param collect_output: Should we collect all the output from the ranks on the master rank?
        :param comm: The MPI communicator used for the parallelization. Default value is None, in which case
            MPI.COMM_WORLD is used

        """
        if not is_mpi_available():
            raise ValueError("MPI is not available. MPI is required for parallel execution.")
        self.task_function = task_function
        self.collect_output = collect_output
        self.schedule = schedule
        self.split_axes = split_axes
        self.main_data = main_data
        self.main_data_param_name = main_data_param_name
        self.task_function_params = task_function_params
        if self.task_function_params is None:
            self.task_function_params = {}
        self.root = root
        self.result = None
        self.comm = get_comm_world() if comm is None else comm

    def run(self):
        """
        Call this function to run the function in parallel
        """
        if self.schedule == self.SCHEDULES['DYNAMIC']:
            return self.run_dynamic()
        elif self.schedule == self.SCHEDULES['STATIC_1D']:
            return self.run_static_1D()

    def run_static_1D(self):
        """
        Run the task function using a static task decomposition schema.

        The data is divided into sub-blocks along the largest split_axis

        :return: i) The result from the local execution of the task_function and
            ii) The selection of the sub-block that was processed. If collect_data
            is True and we are the self.root rank, then the result will be a list
            of the results from all ranks and the sub_block index will be a list of
            all the sub_block indexes used.
        """
        import time
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
        my_block = [slice(None)] * len(self.main_data.shape)
        start_index = rank * block_size
        stop_index = start_index + block_size
        if rank == (size-1):
            if stop_index != split_axis_size:
                stop_index = split_axis_size
        my_block[axes_sort_index[0]] = slice(start_index, stop_index)
        my_block = tuple(my_block)
        print "Rank: " + str(rank) + " Block: " + str(my_block)

        # Execute the task_function on the given data block
        task_params = self.task_function_params
        task_params[self.main_data_param_name] = self.main_data[my_block]
        self.result = self.task_function(**task_params)

        end_time = time.time()
        run_time = end_time - start_time
        print "TIME FOR PROCESSING THE DATA BLOCK: " + str(run_time)

        # Collect the output
        start_time = time.time()
        if self.collect_output:
            collected_data = self.comm.gather(self.result, root=self.root)
            collected_blocks = self.comm.gather(my_block, root=self.root)
            if rank == self.root:
                self.result = collected_data
                my_block = collected_blocks
        end_time = time.time()
        run_time = end_time - start_time
        print "TIME FOR COLLECTING DATA FROM ALL TASKS: " + str(run_time)

        # Return the output
        return self.result, my_block

    def run_dynamic(self):
        """
        Run the task function using dynamic task scheduling.

        The root rank divides the data into sub-tasks and sends the tasks to available MPI
        processes on request.

        :return: The result from all local executions of the tast_function
        """
        import time
        rank = get_rank()
        size = get_size()

        if size < 2:
            warnings.warn('DYNAMIC task scheduling requires at least 2 MPI ranks. Using STATIC scheduling instead.')
            return self.run_static_1D()

        # We are the controlling rank
        if rank == self.root:
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
                if (block_index % 1000) == 0:
                    print (block_index, total_num_subblocks, request_rank)
            end_time = time.time()
            run_time = end_time - start_time
            print "TIME FOR RUNNING TASK FUNCTION: " + str(run_time)
            start_time = time.time()
            print "FINALIZING AND COLLECTING DATA"
            # Terminate all ranks and receive all data from the different ranks if requested
            all_ranks_status = np.zeros(size, 'bool')
            all_ranks_status[self.root] = True
            collected_data = []
            block_selections = []
            while not np.all(all_ranks_status):

                request_rank = self.comm.recv(source=MPI.ANY_SOURCE, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                self.comm.send((None, None), dest=request_rank, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                if self.collect_output:
                    print "COLLECTING DATA FROM: " + str(request_rank)
                    collected_data += self.comm.recv(source=request_rank, tag=self.MPI_MESSAGE_TAGS['COLLECT_MSG'])
                    block_selections += self.comm.recv(source=request_rank, tag=self.MPI_MESSAGE_TAGS['COLLECT_MSG'])
                    print "SIZE: " + str((sys.getsizeof(collected_data), sys.getsizeof(block_selections)))
                all_ranks_status[request_rank] = True
            self.result = collected_data
            end_time = time.time()
            run_time = end_time - start_time
            print "TIME FOR COLLECTING DATA FROM ALL TASKS: " + str(run_time)
            print "DONE COLLECTING DATA"

        # We are a rank that has to run tasks
        else:
            # Request a new data block
            self.result = []
            block_selections = []

            while True:
                self.comm.send(rank, dest=self.root, tag=self.MPI_MESSAGE_TAGS['RANK_MSG'])
                block_index, block_selection = self.comm.recv(source=self.root, tag=self.MPI_MESSAGE_TAGS['BLOCK_MSG'])
                if block_index is None:
                    if self.collect_output:
                        self.comm.send(self.result, dest=self.root, tag=self.MPI_MESSAGE_TAGS['COLLECT_MSG'])
                        self.comm.send(block_selections, dest=self.root, tag=self.MPI_MESSAGE_TAGS['COLLECT_MSG'])
                        # del self.result
                        # self.result=None
                        # del block_selections
                        # block_selections = None
                    break
                # Execute the task_function on the given data block
                task_params = self.task_function_params
                task_params[self.main_data_param_name] = self.main_data[block_selection]
                self.result.append(self.task_function(**task_params))
                block_selections.append(block_selection)

        # Return the result
        return self.result, block_selections



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





