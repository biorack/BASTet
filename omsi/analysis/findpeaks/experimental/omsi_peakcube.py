from time import time, ctime
import sys

import numpy as np

from omsi.analysis.base import analysis_base
from omsi.dataformat.omsi_file.main_file import omsi_file


class omsi_peakcube(analysis_base):
    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""
        super(omsi_peakcube, self).__init__()
        self.analysis_identifier = name_key
        self.data_names = ['npg_peak_cube_mz', 'npg_peak_mz']
        #self.parameters = ['peaksBins', 'peaksIntensities', 'peaksArrayIndex', 'peaksMZdata', 'HCpeaksLabels',
        #                   'HCLabelsList']
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='peakBins',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='peaksIntensities',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='peaksArrayIndex',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='peaksMZdata',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='HCpeaksLabels',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='HCLabelsList',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)

    # ------------------------ viewer functions start ----------------------

    @classmethod
    def v_qslice(cls, analysis_object, z, viewer_option=0):
        """Implement support for qslice URL requests for the viewer"""
        from omsi.shared.data_selection import selection_string_to_object

        if viewer_option == 0:
            dataset = analysis_object['npg_peak_cube_mz']
            try:
                z_select = selection_string_to_object(selection_string=z)
                data = dataset[:, :, z_select]
                return data
            except:
                print "Global peak selection failed"
                return None
        elif viewer_option >= 0:
            return super(omsi_peakcube, cls).v_qslice(analysis_object, z, viewer_option - 1)
        else:
            return None

    @classmethod
    def v_qspectrum(cls, analysis_object, x, y, viewer_option=0):
        """Implement support for qspectrum URL requests for the viewer"""
        # Get the h5py dataset with the peak_cube data
        data = None
        customMZ = None
        if viewer_option == 0:
            from omsi.shared.data_selection import check_selection_string, selection_type, \
                selection_string_to_object

            dataset = analysis_object['npg_peak_cube_mz']
            x_select = selection_string_to_object(selection_string=x)
            y_select = selection_string_to_object(selection_string=y)
            if (check_selection_string(x) == selection_type['indexlist']) and (
                check_selection_string(y) == selection_type['indexlist']):
                #The peak-cube data is usually small enough. To handle the multiple list selection case
                #we here just load the full data cube and use numpy to do the subselection. Note, this
                #version would work for all selection types but we would like to avoid loading the
                #full data if we don't have to.
                data = dataset[:][x_select, y_select, :]
            else:
                data = dataset[x_select, y_select, :]
            #Return the spectra and indicate that no customMZ data values (i.e. None) are needed
            return data, None
        elif viewer_option > 0:
            return super(omsi_peakcube, cls).v_qspectrum(analysis_object, x, y, viewer_option - 1)

        return data, customMZ

    @classmethod
    def v_qmz(cls, analysis_object, qslice_viewer_option=0, qspectrum_viewer_option=0):
        """Implement support for qmz URL requests for the viewer"""
        mzSpectra = None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        peak_cube_shape = analysis_object['npg_peak_cube_mz'].shape
        valuesX = range(0, peak_cube_shape[0])
        labelX = 'pixel index X'
        valuesY = range(0, peak_cube_shape[1])
        labelY = 'pixel index Y'
        if len(peak_cube_shape) > 3:
            valuesZ = range(0, peak_cube_shape[2])
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None
        # We do not need to handle the qslice_viewer_option separately here since there is only one option right now
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  #Loadings
            mzSpectra = analysis_object['npg_peak_mz'][:]
            labelSpectra = "m/z"
            mzSlice = None
            labelSlice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_peakcube, cls).v_qmz(analysis_object, qslice_viewer_option - 1, qspectrum_viewer_option - 1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option > 0:
            mzSpectra = analysis_object['npg_peak_mz'][:]
            labelSpectra = "m/z"
            tempA, tempB, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_peakcube, cls).v_qmz(analysis_object, 0, qspectrum_viewer_option - 1)
        elif qspectrum_viewer_option > 0 and qslice_viewer_option == 0:
            mzSlice = analysis_object['npg_peak_mz'][:]
            labelSlice = "m/z"
            # Ignore the spatial axes. We need to uise the axes of omsi_peakcube
            mzSpectra, labelSpectra, tempA, tempB, vX, lX, vY, lY, vZ, lZ =\
                super(omsi_peakcube, cls).v_qmz(analysis_object, 0, qspectrum_viewer_option - 1)

        return mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_peakcube, cls).v_qspectrum_viewer_options(analysis_object)
        re = ["Peak cube"] + dependent_options
        return re

    @classmethod
    def v_qslice_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_peakcube, cls).v_qslice_viewer_options(analysis_object)
        re = ["Peak cube"] + dependent_options
        return re

    # ------------------------ viewer functions end ------------------------

    # main PeakCube function
    # omsi_peakcube_exec(self, peaksBins, peaksIntensities, peaksArrayIndex, peaksMZdata, HCpeaksLabels, HCLabelsList, npgdata_dependency = None, lpfdata_dependency = None):
    def execute_analysis(self):
        print "Generating the peak cube..."

        #Copy parameters to local variables for convenience
        peaksBins = self['peaksBins']
        peaksIntensities = self['peaksIntensities']
        peaksArrayIndex = self['peaksArrayIndex']
        peaksMZdata = self['peaksMZdata']
        HCpeaksLabels = self['HCpeaksLabels']
        HCLabelsList = self['HCLabelsList']

        myPeakCube = self.getPeakCube(peaksIntensities, peaksArrayIndex, HCpeaksLabels, HCLabelsList)
        print "\ndone!"
        print "Calculating global m\z..."
        myGlobalMz = self.getGlobalMz(peaksBins, peaksMZdata, HCpeaksLabels, HCLabelsList)
        print "\ndone!"

        print "Collecting data into HDF5..."

        #Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
        #We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        #handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        #to ensure a consitent behavior we convert the values directly here

        #Clear any previously stored analysis data
        self.clear_analysis_data()
        self.clear_parameter_data()

        # Collect peak cube into hdf5
        self['npg_peak_cube_mz'] = np.asarray(myPeakCube)
        self['npg_peak_mz'] = np.asarray(myGlobalMz)

        #		#Save the analysis dependencies to the __dependency_list so that the data can be saved automatically by the omsi HDF5 file API
        #		if npgdata_dependency is not None :
        #			# if isinstance( npgdata_dependency , dependency_dict ) :
        #			# 	#Add the depency as given by the user
        #			# 	self.add_dependency_data( npgdata_dependency )
        #			# elif isinstance( lpfdata_dependency, dependency_dict):
        #			# 	self.add_dependency_data( lpfdata_dependency )
        #
        #			# else :
        #			# 	#The user only gave us the object that we depend on so we need to construct the
        #			self.add_dependency_data( dependency_dict( param_name = 'HCpeaksLabels', link_name='HCpeaksLabels', omsi_object=npgdata_dependency['npghc_peaks_labels'], selection=None ) )
        #			self.add_dependency_data( dependency_dict( param_name = 'HCLabelsList', link_name='HCLabelsList', omsi_object=npgdata_dependency['npghc_labels_list'], selection=None ) )
        #
        #		if lpfdata_dependency is not None :
        #
        #			self.add_dependency_data( dependency_dict( param_name = 'peaksBins', link_name='peaksBins', omsi_object=lpfdata_dependency['LPF_Peaks_MZ'], selection=None ) )
        #			self.add_dependency_data( dependency_dict( param_name = 'peaksIntensities', link_name='peaksIntensities', omsi_object=lpfdata_dependency['LPF_Peaks_Vals'], selection=None ) )
        #			self.add_dependency_data( dependency_dict( param_name = 'peaksArrayIndex', link_name='peaksArrayIndex', omsi_object=lpfdata_dependency['LPF_Peaks_ArrayIndex'], selection=None ) )

        print "Collecting done."
        print "--- finished ---"
        return self['npg_peak_cube_mz'], self['npg_peak_mz']

    def record_execute_analysis_outputs(self, analysis_output):
        """We are recording our outputs manually as part of the execute function"""
        pass

    # +++++++++++++++++++++++ helper functions ++++++++++++++++++++++++++

    def getPeakCube(self, peaksIntensities, peaksArrayIndex, HCpeaksLabels, HCLabelsList):

        pAILast = peaksArrayIndex.shape[0] - 1
        Nx = peaksArrayIndex[pAILast][0] + 1
        Ny = peaksArrayIndex[pAILast][1] + 1
        Nz = len(HCLabelsList)
        PC = np.zeros((Nx * Ny, Nz))

        percentcheck = None
        timekeep1 = time()
        for iNxy in xrange(Nx * Ny):
            percent = int(100. * float(iNxy + 1) / float(Nx * Ny))
            if (percent != percentcheck):
                timer = int(time() - timekeep1)
                print "[", percent, "% -", timer, "s ]\r",
                sys.stdout.flush()
                percentcheck = percent

            pStart = peaksArrayIndex[iNxy][2]

            # look at next arrayindex to see current end,
            # unless its the last arrayindex
            if (iNxy == ((Nx * Ny) - 1)):
                pEnd = HCpeaksLabels.shape[0]
            else:
                pEnd = peaksArrayIndex[iNxy + 1][2]

            if (pStart != pEnd):
                currentInts = peaksIntensities[pStart:pEnd]
                currentLbls = HCpeaksLabels[pStart:pEnd]

                # method 3
                # Does direct index replacement if labels for current coordinate are unique
                uniquelbls = np.unique(currentLbls)
                if (len(uniquelbls) == len(currentLbls)):
                    currentIdxs = currentLbls.astype('int32') - 1
                    PC[iNxy][currentIdxs] = currentInts
                # if not unique then use the max intensity for the same labels
                else:
                    sidxs = np.argsort(currentLbls)
                    srtLbls = currentLbls[sidxs]
                    srtInts = currentInts[sidxs]

                    prevLbl = None
                    for iLbl in xrange(len(currentLbls)):
                        thisLbl = srtLbls[iLbl]
                        thisInt = srtInts[iLbl]
                        myIdx = int(thisLbl - 1)
                        if (thisLbl == prevLbl):
                            prevInt = PC[iNxy][myIdx]
                            myInt = max(prevInt, thisInt)
                        else:
                            myInt = thisInt

                        prevLbl = thisLbl
                        PC[iNxy][myIdx] = myInt

        PC = np.reshape(PC, (Nx, Ny, Nz))
        return PC

    # calculate the global m\z
    def getGlobalMz(self, peaksBins, peaksMZdata, HCpeaksLabels, HCLabelsList):

        import bisect

        timekeep = time()

        labelsArgs = np.argsort(HCpeaksLabels)
        SortedPL = HCpeaksLabels[labelsArgs]
        globalMZ = np.empty_like(HCLabelsList, dtype=float)
        idxLeft = 0

        percentcheck = None
        for i in xrange(len(HCLabelsList)):
            currentlabel = HCLabelsList[i]

            idxRight = bisect.bisect_right(SortedPL, currentlabel)
            currentMZs = peaksMZdata[peaksBins[labelsArgs[idxLeft:idxRight]]]
            globalMZ[i] = np.median(currentMZs)
            idxLeft = idxRight

            percent = int(100. * float(i + 1) / float(len(HCLabelsList)))
            timer = str(int(time() - timekeep))
            if (percent != percentcheck):
                sys.stdout.write("[ " + str(i + 1) + " of " + str(len(HCLabelsList)) + " - " + str(
                    percent) + "% - " + timer + "s ]" + "\r")
                sys.stdout.flush()
                percentcheck = percent

            # end of for i in xrange loop ---

        return globalMZ


def main(argv=None):
    if argv is None:
        argv = sys.argv

        # Check for correct usage
    if len(argv) != 6:
        print "USAGE: Call \"omsi_peakcube OMSI_FILE [expIndex dataIndex LPFanalysisIndex NPGanalysisIndex]   \" "
        print "Generates PeakCube for NPG global peaks"
        exit(0)

    #Read the input arguments
    omsiInFile = argv[1]
    expIndex = int(argv[2])
    dataIndex = int(argv[3])
    LPFanalysisIndex = int(argv[4])
    NPGanalysisIndex = int(argv[5])

    #Open the input HDF5 file
    try:
        omsiFile = omsi_file(omsiInFile)
    except:
        print "Error opening input file:", sys.exc_info()[0]
        exit(0)

    print "Input file: ", omsiInFile
    print "LPF Analysis Index: ", LPFanalysisIndex
    print "NPG Analysis Index: ", NPGanalysisIndex

    #Get the experiment and data
    exp = omsiFile.get_experiment(expIndex)
    data = exp.get_msidata(dataIndex)
    peaksMZdata = data.mz[:]

    LPFanalysis = exp.get_analysis(LPFanalysisIndex)
    NPGanalysis = exp.get_analysis(NPGanalysisIndex)

    print "\n[loading data...]"
    peaksBins = LPFanalysis['LPF_Peaks_MZ'][:]
    peaksIntensities = LPFanalysis['LPF_Peaks_Vals'][:]
    peaksArrayIndex = LPFanalysis['LPF_Peaks_ArrayIndex'][:]
    NPGPL = NPGanalysis['npghc_peaks_labels'][:]
    NPGLL = NPGanalysis['npghc_labels_list'][:]

    print "[done!] lpf data shapes:"
    print "peaksBins shape: ", peaksBins.shape
    print "peaksIntensities shape: ", peaksIntensities.shape
    print "peaksArrayIndex shape: ", peaksArrayIndex.shape
    print "peaksMZdata shape: ", peaksMZdata.shape
    print "NPGPL shape: ", NPGPL.shape
    print "NPGLL shape: ", NPGLL.shape
    sys.stdout.flush()

    # pc
    myPC = omsi_peakcube(name_key="omsi_peakcube_" + str(ctime()))
    print "--- Creating Peak Cube ---"
    #myPC.omsi_peakcube_exec(peaksBins, peaksIntensities, peaksArrayIndex, peaksMZdata, NPGPL, NPGLL)
    myPC.execute(peaksBins=peaksBins,
                 peaksIntensities=peaksIntensities,
                 peaksArrayIndex=peaksArrayIndex,
                 peaksMZdata=peaksMZdata,
                 HCpeaksLabels=NPGPL,
                 HCLabelsList=NPGLL)

    PCm = myPC['npg_peak_cube_mz']
    print "NPG Peak Cube Mzs: \n", PCm.shape, "\n", PCm
    PMz = myPC['npg_peak_mz']
    print "NPG Peak Mz: \n", PMz.shape, "\n", PMz

    print "\nsaving HDF5 analysis..."
    PCanalysis, PCanalysisindex = exp.create_analysis(myPC)
    print "done!"
    print "--- omsi_peakcube complete ---\n"
    print "Peak Cube analysis index:", PCanalysisindex
    print "omsi_peakcube complete for input file: ", omsiInFile, "\n"


# stop python, used for debugging
def stop():
    raw_input("Stop!")


if __name__ == "__main__":
    main()
