from time import time, ctime
import math
import sys
from sys import exit

import numpy as np

from omsi.analysis.base import analysis_base
from omsi.dataformat.omsi_file.main_file import omsi_file

PRECISION = "double"
SLWMIN_ARR_SIZE = 10000  # size of the dequeue for slmin
PKSDET_ARR_SIZE = 4000  # limit of total of peaks that can be stored (C array sizes) per x/y spectra
SET_MEM_CAP = 1200000000  # 1200000000 bytes = 1.12 gb :: max size of data per cl gate
openCLfpFile = "opencl_findpeaks.cl"  # file name of the findpeaks opencl code

# class omsi_findpeaks_local(analysis_base) :
class omsi_lpf(analysis_base):
    def __init__(self, name_key="undefined"):
        super(omsi_lpf, self).__init__()
        self.analysis_identifier = name_key
        # self.parameters = [ 'msidata' , 'mzdata', 'integration_width', 'peakheight', 'slwindow', 'smoothwidth']
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='msidata',
                           help='The MSI dataset to be analyzed',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='mzdata',
                           help='The m/z values for the spectra of the MSI dataset',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='integration_width',
                           help='The window over which peaks should be integrated',
                           dtype=float,
                           default=10,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='peakheight',
                           help='Peak height parameter',
                           dtype=int,
                           default=2,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='slwindow',
                           help='Sliding window parameter',
                           dtype=int,
                           default=100,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='smoothwidth',
                           help='Smooth width parameter',
                           dtype=int,
                           default=3,
                           group=groups['settings'],
                           required=True)
        self.data_names = ['LPF_Peaks_MZ', 'LPF_Peaks_Vals', 'LPF_Peaks_ArrayIndex', 'LPF_indata_mz']


# ------------------------ viewer functions start ----------------------

@classmethod
def v_qslice(cls, analysis_object, z, viewer_option=0):
    """Implement support for qslice URL requests for the viewer"""
    # Use the dependency data for slicing here. We do not have a native option to reconstruct images from local peak finding data
    return super(omsi_lpf, cls).v_qslice(analysis_object, z, viewer_option)


@classmethod
def v_qspectrum(cls, analysis_object, x, y, viewer_option=0):
    """Implement support for qspectrum URL requests for the viewer"""
    # Retrieve the h5py objects for the requried datasets from the local peak finding
    if viewer_option == 0:

        from omsi.shared.data_selection import check_selection_string, selection_type, selection_to_indexlist
        import numpy as np

        peak_mz = analysis_object['LPF_Peaks_MZ']
        peak_values = analysis_object['LPF_Peaks_Vals']
        arrayIndices = analysis_object['LPF_Peaks_ArrayIndex'][:]
        indata_mz = analysis_object['LPF_indata_mz']
        #Determine the shape of the original raw data
        if (indata_mz is None) or (arrayIndices is None):
            return None, None
        Nx = arrayIndices[:, 0].max()
        Ny = arrayIndices[:, 1].max()
        NMZ = indata_mz.shape[0]
        numSpectra = arrayIndices.shape[0]
        #Determine the size of the selection and the set of selected items
        xList = selection_to_indexlist(x, Nx)
        yList = selection_to_indexlist(y, Ny)
        if (check_selection_string(x) == selection_type['indexlist']) and (
            check_selection_string(y) == selection_type['indexlist']):
            if len(xList) == len(yList):
                items = [(xList[i], yList[i]) for i in xrange(0, len(xList))]
            else:
                return None, None
        else:
            items = [0] * (len(xList) * len(yList))
            index = 0
            for xi in xList:
                for yi in yList:
                    items[index] = (xi, yi )
                    index = index + 1

        shapeX = len(items)
        shapeY = 1
        shapeZ = NMZ
        #Initalize the data cube to be returned
        data = np.zeros((shapeX, shapeY, shapeZ), dtype=peak_values.dtype)
        #Fill the non-zero locations for the data cube with data
        for ni in xrange(0, len(items)):
            currentIndex = (items[ni][0] * Ny + items[ni][1])
            currentDX = ni
            currentDY = 0
            startIndex = arrayIndices[currentIndex][2]
            if currentIndex < numSpectra:
                endIndex = arrayIndices[(currentIndex + 1)][2]
            else:
                endIndex = peak_values.size
            if startIndex != endIndex:
                tempValues = peak_values[startIndex: endIndex]
                tempMZ = peak_mz[startIndex: endIndex]
                data[currentDX, currentDY, tempMZ] = tempValues
            else:
                #The start and end index may be the same in case that no peaks for found for the given spectrum
                #The data is already initalized to 0 so there is nothing to do here
                pass

        if len(items) == 1:
            data = data.reshape((shapeX, shapeZ))

        #Return the spectra and indicate that no customMZ data values (i.e. None) are needed
        return data, None

    elif viewer_option > 0:
        return super(omsi_lpf, cls).v_qspectrum(analysis_object, x, y, viewer_option - 1)
    else:
        return None, None


@classmethod
def v_qmz(cls, analysis_object, qslice_viewer_option=0, qspectrum_viewer_option=0):
    """Implement support for qmz URL requests for the viewer"""
    mzSpectra = None
    labelSpectra = None
    mzSlice = None
    labelSlice = None
    array_indices = analysis_object['peak_arrayindex'][:]
    x_size = array_indices[:, 0].max()+1
    y_size = array_indices[:, 1].max()+1
    valuesX = range(0, x_size)
    labelX = 'pixel index X'
    valuesY = range(0, y_size)
    labelY = 'pixel index Y'
    if array_indices.shape[1] > 2:
        z_size = array_indices[:, 2].max()+1
        valuesZ = range(0, z_size)
        labelZ = 'pixel index Z'
    else:
        valuesZ = None
        labelZ = None

    # We do not have native option for qslice, so we rely on the input data in all cases
    if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  #Loadings
        mzSpectra = analysis_object['LPF_indata_mz'][:]
        labelSpectra = "m/z"
        mzSlice = None
        labelSlice = None
    elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
        mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ  = \
            super(omsi_lpf, cls).v_qmz(analysis_object, qslice_viewer_option, qspectrum_viewer_option - 1)
    elif qspectrum_viewer_option == 0 and qslice_viewer_option >= 0:
        mzSpectra = analysis_object['LPF_indata_mz'][:]
        labelSpectra = "m/z"
        tempA, tempB, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ  = \
            super(omsi_lpf, cls).v_qmz(analysis_object, 0, qspectrum_viewer_option)

    return mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ


@classmethod
def v_qspectrum_viewer_options(cls, analysis_object):
    """Define which viewer_options are supported for qspectrum URL's"""
    dependent_options = super(omsi_lpf, cls).v_qspectrum_viewer_options(analysis_object)
    re = ["Local peaks"] + dependent_options
    return re


@classmethod
def v_qslice_viewer_options(cls, analysis_object):
    """Define which viewer_options are supported for qspectrum URL's"""
    return super(omsi_lpf, cls).v_qslice_viewer_options(analysis_object)


# ------------------------ viewer functions end ------------------------
def execute_analysis(self):
    import pyopencl as cl

    #Copy parameters to local variables for convenience
    msidata = self['msidata']
    mzdata = self['mzdata']
    integration_width = self['integration_width']
    peakheight = self['peakheight']
    slwindow = self['slwindow']
    smoothwidth = self['smoothwidth']

    #Determine the data dimensions
    Nx = msidata.shape[0]
    Ny = msidata.shape[1]
    Nz = msidata.shape[2]
    print "Data shape: ", msidata.shape
    print "[LPF]"
    print "LPF Parameters:"
    print "peakheight =", peakheight, "slwindow =", slwindow, "smoothwidth =", smoothwidth
    sys.stdout.flush()

    timebefore = time()  # timekeeping

    peak_values = []
    peak_MZ = []
    peak_arrayindex = np.zeros(shape=(Nx * Ny, 3), dtype='int64')
    pkcount = 0
    cindex = 0
    pcount = 0

    # float64 = double, float32 = float
    if (PRECISION == "double"):
        NP_TYPE = np.float64
        mytypedef = "double"
    else:
        NP_TYPE = np.float32
        mytypedef = "float"

    cl_info = cl.get_platforms()[0]
    cl_info = cl_info.get_devices(device_type=cl.device_type.GPU)
    cl_info = len(cl_info)

    # memory cap based on GPUs ammount
    print "# of GPUs: ", cl_info
    MEM_CAP = SET_MEM_CAP * cl_info * .90

    # x/y split due to gpu memory constraint
    xsplit = 1
    ysplit = Ny
    for Nxi in xrange(1, Nx + 1):
        splitarr = np.empty([Nxi, Ny, Nz], dtype=NP_TYPE)
        if (splitarr.nbytes < MEM_CAP):
            xsplit = Nxi
        else:
            break

    if (xsplit == 1):
        ysplit = 1
        for Nyi in xrange(1, Ny + 1):
            splitarr = np.empty([1, Nyi, Nz], dtype=NP_TYPE)
            if (splitarr.nbytes < MEM_CAP):
                ysplit = Nyi
            else:
                break
        for Nyi in xrange(1, Ny + 1):
            ydiv = Ny / float(Nyi)
            ydiv = int(math.ceil(ydiv))
            if (ysplit >= ydiv):
                ysplit = ydiv
                break

    print "xsplit: ", xsplit, " ysplit: ", ysplit

    xstart = ystart = 0
    xend = xamt = xsplit
    yend = yamt = ysplit

    xcount = Nx / float(xsplit)
    xcount = int(math.ceil(xcount))
    ycount = Ny / float(ysplit)
    ycount = int(math.ceil(ycount))

    for xc in xrange(xcount):
        for yc in xrange(ycount):
            print "#: ", pcount + 1, "/", xcount * ycount
            sys.stdout.flush()

            # opencl gate
            [cl_maxtabVal, cl_maxtabPos, cl_maxtabTot] = self.cl_peakfind(msidata[xstart:xend, ystart:yend, :],
                                                                          smoothwidth, slwindow, peakheight)

            # building peak structures
            for x in xrange(xamt):
                for y in xrange(yamt):

                    totindex = y + ((yend - ystart) * x)
                    for t in xrange(cl_maxtabTot[totindex]):
                        peak_values.append(cl_maxtabVal[x, y, t])
                        peak_MZ.append(cl_maxtabPos[x, y, t])

                    peak_arrayindex[pkcount, 0] = x + (xsplit * xc)
                    peak_arrayindex[pkcount, 1] = y + (ysplit * yc)
                    peak_arrayindex[pkcount, 2] = cindex

                    pkcount += 1
                    cindex += cl_maxtabTot[totindex]

            ystart = yend
            yend += ysplit
            if (yend > Ny):
                yend = Ny
                yamt = yend - ystart
            if (ystart >= Ny):
                ystart = 0
                yend = ysplit
                yamt = ysplit
                xstart = xend
                xend += xsplit

            if (xend > Nx):
                xend = Nx
                xamt = xend - xstart
            if (xstart >= xend):
                break

            pcount += 1

    print "peak_arrayindex shape: ", np.array(peak_arrayindex).shape
    print "peak_values: ", np.array(peak_values)
    print "peak_values shape: ", np.array(peak_values).shape
    print "peak_MZ: ", np.array(peak_MZ)
    print "peak_MZ shape: ", np.array(peak_MZ).shape

    timetotal = time() - timebefore  # timekeeping
    print "time (s): ", timetotal  # timekeeping

    print "Collecting data into HDF5..."

    #Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
    #We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
    #handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
    #to ensure a consitent behavior we convert the values directly here

    #Save the analysis data to the __data_list so that the data can be saved automatically by the omsi HDF5 file API
    return np.asarray(peak_MZ), np.asarray(peak_values), np.asarray(peak_arrayindex), mzdata[:]
    # self['LPF_Peaks_MZ'] = np.asarray( peak_MZ )
    # self['LPF_Peaks_Vals'] = np.asarray( peak_values )
    # self['LPF_Peaks_ArrayIndex'] = np.asarray( np.asarray( peak_arrayindex ) )
    # self['LPF_Peaks_Vals'] = np.asarray( peak_values )
    # self['LPF_indata_mz'] = mzdata[:]

    print "Collecting done."
    print "--- finished ---"


def cl_peakfind(self, msidt, smoothsize, slwindow, peakheight):
    import pyopencl as cl

    # float64 = double, float32 = float
    if (PRECISION == "double"):
        NP_TYPE = np.float64
        mytypedef = "double"
    else:
        NP_TYPE = np.float32
        mytypedef = "float"

    dev_amt = len(cl.get_platforms()[0].get_devices(device_type=cl.device_type.GPU))

    msidt2d = np.reshape(msidt, (msidt.shape[0] * msidt.shape[1], msidt.shape[2]))
    total_split = int(math.ceil(msidt2d.shape[0] / float(dev_amt)))
    current_split = 0

    in_data = []
    xsize = []
    ysize = []

    for d in xrange(dev_amt):
        in_data.append(np.array(msidt2d[current_split:total_split, :], dtype=NP_TYPE))
        xsize.append(np.int32(in_data[d].shape[0]))
        ysize.append(np.int32(in_data[d].shape[1]))
        print "sizes", d + 1, ":", xsize[d], ysize[d], "  ",
        current_split = total_split
        total_split += total_split
        if (total_split > msidt2d.shape[0]):
            total_split = msidt2d.shape[0]

    # ------------- opencl init [s] -------------
    dstr = "// opencl code \n"
    if (mytypedef == "double"):
        dstr += "#pragma OPENCL EXTENSION cl_khr_fp64: enable \n"
    dstr += "typedef " + mytypedef + " clfloat; \n"

    mf = cl.mem_flags
    f = open(openCLfpFile, 'r')
    fstr = "".join(f.readlines())
    fstr = dstr + fstr
    f.close()

    cl_pl = cl.get_platforms()[0]
    cl_dev = cl_pl.get_devices(device_type=cl.device_type.GPU)

    ctx = []
    queue = []
    prg = []
    for d in xrange(dev_amt):
        ctx.append(cl.Context(devices=[cl_dev[d]]))
        queue.append(cl.CommandQueue(ctx[d]))
        prg.append(cl.Program(ctx[d], fstr).build())

    # ------------- opencl init [e] -------------

    # ------ smoothgaussian [s] ------
    x = np.array(range(-3 * smoothsize, 3 * smoothsize))
    g = np.exp(-(x ** 2) / (2.0 * float(smoothsize) ** 2))
    g = g / g.sum()
    g = np.array(g, dtype=NP_TYPE)

    smoothsize = np.int32(smoothsize)

    input_buf = []
    g_buf = []
    output_buf = []
    global_size = []
    local_size = None
    for d in xrange(dev_amt):
        input_buf.append(cl.Buffer(ctx[d], mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=in_data[d]))
        g_buf.append(cl.Buffer(ctx[d], mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=g))
        output_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, in_data[d].nbytes))
        global_size.append((in_data[d].shape[0],))
        kernelargs = (    input_buf[d],
                          g_buf[d],
                          smoothsize,
                          xsize[d], ysize[d],
                          output_buf[d])
        prg[d].clSmoothed(queue[d], global_size[d], local_size, *(kernelargs))

    for d in xrange(dev_amt):
        queue[d].finish()

        cl.enqueue_copy(queue[d], in_data[d], output_buf[d])

        input_buf[d].release()
        g_buf[d].release()
        output_buf[d].release()

    # ------ smoothgaussian [e] ------
    print "\nsmooth - done // ",
    sys.stdout.flush()
    # ------ sliding window [s] ------
    slwindow = np.int32(slwindow)
    window_arr_size = np.int32(SLWMIN_ARR_SIZE)

    input_buf = []
    value_buf = []
    index_buf = []
    output_buf = []

    for d in xrange(dev_amt):
        dummyVal = np.empty_like(in_data[d][:, 0:SLWMIN_ARR_SIZE], dtype=NP_TYPE)
        dummyPos = np.empty_like(dummyVal, dtype=np.int32)

        input_buf.append(cl.Buffer(ctx[d], mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=in_data[d]))
        value_buf.append(cl.Buffer(ctx[d], mf.READ_WRITE, dummyVal.nbytes))
        index_buf.append(cl.Buffer(ctx[d], mf.READ_WRITE, dummyPos.nbytes))
        output_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, in_data[d].nbytes))

        kernelargs = (    input_buf[d],
                          value_buf[d], index_buf[d],
                          window_arr_size,
                          slwindow,
                          xsize[d], ysize[d],
                          output_buf[d])

        prg[d].clSliding_window_minimum(queue[d], global_size[d], local_size, *(kernelargs))

    for d in xrange(dev_amt):
        queue[d].finish()

        cl_slmin = np.empty_like(in_data[d])
        cl.enqueue_copy(queue[d], cl_slmin, output_buf[d])

        in_data[d] = in_data[d] - cl_slmin

        input_buf[d].release()
        value_buf[d].release()
        index_buf[d].release()
        output_buf[d].release()

    # ------ sliding window [e] ------
    print "sliding - done // ",
    sys.stdout.flush()
    # --------- peakdet [s] ----------

    my_size_limit = PKSDET_ARR_SIZE  # limit of total of peaks that can be stored (C array sizes)
    cl_maxtabTotal = np.int32(my_size_limit)
    cl_mintabTotal = np.int32(my_size_limit)
    peakheight = np.int32(peakheight)

    input_buf = []
    maxtabVal_buf = []
    maxtabPos_buf = []
    maxtabTot_buf = []
    mintabVal_buf = []
    mintabPos_buf = []
    mintabTot_buf = []
    cl_maxtabVal = []
    cl_maxtabPos = []
    cl_maxtabTot = []
    cl_mintabVal = []
    cl_mintabPos = []
    cl_mintabTot = []

    for d in xrange(dev_amt):
        dummyVal = np.empty_like(in_data[d][:, 0:my_size_limit], dtype=NP_TYPE)
        dummyPos = np.empty_like(dummyVal, dtype=np.int32)
        dummyTot = np.array(dummyVal[:, 0], dtype=np.int32)

        input_buf.append(cl.Buffer(ctx[d], mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=in_data[d]))
        maxtabVal_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyVal.nbytes))
        maxtabPos_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyPos.nbytes))
        maxtabTot_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyTot.nbytes))
        mintabVal_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyVal.nbytes))
        mintabPos_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyPos.nbytes))
        mintabTot_buf.append(cl.Buffer(ctx[d], mf.WRITE_ONLY, dummyTot.nbytes))

        kernelargs = (    input_buf[d],
                          maxtabVal_buf[d], maxtabPos_buf[d], maxtabTot_buf[d],
                          mintabVal_buf[d], mintabPos_buf[d], mintabTot_buf[d],
                          peakheight,
                          xsize[d], ysize[d],
                          cl_maxtabTotal, cl_mintabTotal    )

        prg[d].clPeakdet(queue[d], global_size[d], local_size, *(kernelargs))

        cl_maxtabVal.append(np.empty_like(dummyVal))
        cl_maxtabPos.append(np.empty_like(dummyPos))
        cl_maxtabTot.append(np.empty_like(dummyTot))
        cl_mintabVal.append(np.empty_like(dummyVal))
        cl_mintabPos.append(np.empty_like(dummyPos))
        cl_mintabTot.append(np.empty_like(dummyTot))

    for d in xrange(dev_amt):
        queue[d].finish()

        cl.enqueue_copy(queue[d], cl_maxtabVal[d], maxtabVal_buf[d])
        cl.enqueue_copy(queue[d], cl_maxtabPos[d], maxtabPos_buf[d])
        cl.enqueue_copy(queue[d], cl_maxtabTot[d], maxtabTot_buf[d])
        cl.enqueue_copy(queue[d], cl_mintabVal[d], mintabVal_buf[d])
        cl.enqueue_copy(queue[d], cl_mintabPos[d], mintabPos_buf[d])
        cl.enqueue_copy(queue[d], cl_mintabTot[d], mintabTot_buf[d])

        input_buf[d].release()
        maxtabVal_buf[d].release()
        maxtabPos_buf[d].release()
        maxtabTot_buf[d].release()
        mintabVal_buf[d].release()
        mintabPos_buf[d].release()
        mintabTot_buf[d].release()

    # --------- peakdet [e] ----------
    print "peakdet - done"
    sys.stdout.flush()

    cl_maxtabVal = np.concatenate(cl_maxtabVal[:])
    cl_maxtabPos = np.concatenate(cl_maxtabPos[:])
    cl_maxtabTot = np.concatenate(cl_maxtabTot[:])

    cl_maxtabVal = np.reshape(cl_maxtabVal, (msidt.shape[0], msidt.shape[1], cl_maxtabVal[1]))
    cl_maxtabPos = np.reshape(cl_maxtabPos, (msidt.shape[0], msidt.shape[1], cl_maxtabPos[1]))
    # cl_maxtabPos = np.reshape(cl_maxtabPos,(msidt.shape[0], msidt.shape[1]))

    return cl_maxtabVal, cl_maxtabPos, cl_maxtabTot


def main(argv=None):
    """Then main function"""

    if argv is None:
        argv = sys.argv

        # Check for correct usage
    if len(argv) < 2:
        print "USAGE: Call \"omsi_lpf OMSI_FILE [expIndex dataIndex peakHeight]	\" "
        print "Local peakfinding using opencl"
        exit(0)

    #Read the input arguments
    omsiInFile = argv[1]

    sys.stdout.flush()
    expIndex = 0
    dataIndex = 0
    mypeakheight = 10

    if len(argv) == 4:
        expIndex = int(argv[2])
        dataIndex = int(argv[3])

    if len(argv) == 5:
        expIndex = int(argv[2])
        dataIndex = int(argv[3])
        mypeakheight = int(argv[4])

    #Open the input HDF5 file
    try:
        omsiFile = omsi_file(omsiInFile, 'a')
    except:
        print "Unexpected openeing the input file:", sys.exc_info()[0]
        exit(0)

    print "Input file: ", omsiInFile

    #Get the experiment and dataset
    exp = omsiFile.get_experiment(expIndex)
    data = exp.get_msidata(dataIndex)
    mzdata = data.mz[:]

    #Execute the peak finding
    myLPF = omsi_lpf(name_key="omsi_lpf_" + str(ctime()))
    print "--- Executing LPF ---"
    myLPF.execute(data, mzdata, peakheight=mypeakheight)
    print "\n\nGetting peak finding analysis results"
    pmz = myLPF['LPF_Peaks_MZ']
    print "pmz:\n", pmz
    pv = myLPF['LPF_Peaks_Vals']
    print "pv:\n", pv
    pai = myLPF['LPF_Peaks_ArrayIndex']
    print "pai:\n", pai

    print "\nsaving HDF5 analysis..."
    analysis, analysisindex = exp.create_analysis(myLPF)
    print "done!"

    print "lpf analysis index:", analysisindex
    print "omsi_lpf complete for input file: ", omsiInFile, "\n"


if __name__ == "__main__":
    main()

