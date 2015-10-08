from time import time, ctime
import sys
from sys import exit

import numpy as np

from omsi.analysis.base import analysis_base
from omsi.dataformat.omsi_file.main_file import omsi_file


# ----- hier clustering imports
# import fastcluster
# import scipy
# import scipy.cluster.hierarchy as hier

# Node class for using union/find functions on LabelsLists Node elements
# import itertools
class Node:
    def __init__(self, label):
        self.label = label

    def __str__(self):
        return str(self.label)


class omsi_npg(analysis_base):
    # class attributes: LPF data
    # set on start of omsi_npg_exec
    ##These have been replace by the parameters infrastrucutre
    # peaksBins = None
    #peaksIntensities = None
    #peaksArrayIndex = None
    #peaksMZdata = None
    #peaksMZ = None

    def __init__(self, name_key="undefined"):
        """Initalize the basic data members"""
        super(omsi_npg, self).__init__()
        self.analysis_identifier = name_key
        self.data_names = ['npg_labels_medianmz', 'npg_labels_list', 'npg_peaks_labels', 'npghc_tc_labels',
                           'npghc_labels_list', 'npghc_peaks_labels']
        # self.parameters = ['peaksBins', 'npg_peaks_Intensities', 'npg_peaks_ArrayIndex', 'peaksMZdata', 'peaksMZ',
        #                   'npg_cluster_method', 'npg_split_max', 'npg_fasterhc_flag', 'npg_mz_threshold',
        #                   'npg_cluster_treecut']
        dtypes = self.get_default_dtypes()
        groups = self.get_default_parameter_groups()
        self.add_parameter(name='peakBins',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='npg_peaks_Intensities',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='npg_peaks_ArrayIndex',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='peaksMZdata',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='peaksMZ',
                           help='',
                           dtype=dtypes['ndarray'],
                           group=groups['input'],
                           required=True)
        self.add_parameter(name='npg_cluster_method',
                           help='Clustering method to be used',
                           dtype=unicode,
                           default='median',
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='npg_split_max',
                           help='',
                           dtype=int,
                           default=1000,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='npg_fasterhc_flag',
                           help='',
                           dtype=int,
                           default=3,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='npg_mz_threshold',
                           help='',
                           dtype=float,
                           default=0.05,
                           group=groups['settings'],
                           required=True)
        self.add_parameter(name='npg_cluster_treecut',
                           help='Threshold where the hierarichal clustering tree should be cur',
                           dtype=float,
                           default=0.1,
                           group=groups['settings'],
                           required=True)

    #
    # ------------------------ viewer functions start ----------------------

    @classmethod
    def v_qslice(cls, analysis_object, z, viewer_option=0):
        """Implement support for qslice URL requests for the viewer"""
        from omsi.shared.data_selection import selection_string_to_object

        if viewer_option == 0:
            pl = analysis_object['npg_peaks_labels']
            ll = analysis_object['npg_labels_list']
            pai = analysis_object['npg_peaks_ArrayIndex']
            pints = analysis_object['npg_peaks_Intensities']
            z_select = selection_string_to_object(selection_string=z)
            dataset = cls.getnpgimage(pl, ll, pai, pints, z_select)
            try:
                return dataset[:, :]
            except:
                print "NPG peak selection failed"
                return None
        elif viewer_option >= 0:
            return super(omsi_npg, cls).v_qslice(analysis_object, z, viewer_option - 1)
        else:
            return None

    @classmethod
    def v_qspectrum(cls, analysis_object, x, y, viewer_option=0):
        """Implement support for qspectrum URL requests for the viewer"""
        #Get the h5py dataset with the peak_cube data
        data = None
        customMZ = None
        if viewer_option == 0:
            from omsi.shared.data_selection import selection_string_to_object

            pl = analysis_object['npg_peaks_labels']
            ll = analysis_object['npg_labels_list']
            pai = analysis_object['npg_peaks_ArrayIndex']
            pints = analysis_object['npg_peaks_Intensities']
            x_select = selection_string_to_object(selection_string=x)
            y_select = selection_string_to_object(selection_string=y)
            dataset = cls.getnpgspec(pl, ll, pai, pints, x_select, y_select)

            try:
                return dataset[:]
            except:
                print "NPG spectrum selection failed"
                return None

            #Return the spectra and indicate that no customMZ data values (i.e. None) are needed
            return data, None
        elif viewer_option > 0:
            return super(omsi_npg, cls).v_qspectrum(analysis_object, x, y, viewer_option - 1)

        return data, customMZ

    @classmethod
    def v_qmz(cls, analysis_object, qslice_viewer_option=0, qspectrum_viewer_option=0):
        """Implement support for qmz URL requests for the viewer"""
        mzSpectra = None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        peaksArrayIndex = analysis_object['npg_peaks_ArrayIndex']
        pAILast = peaksArrayIndex.shape[0] - 1
        x_size = peaksArrayIndex[pAILast][0] + 1
        y_size = peaksArrayIndex[pAILast][1] + 1
        valuesX = range(0, x_size)
        labelX = 'pixel index X'
        valuesY = range(0, y_size)
        labelY = 'pixel index Y'
        if peaksArrayIndex.shape[1] > 2:
            z_size = peaksArrayIndex[pAILast][2] + 1
            valuesZ = range(0, z_size)
            labelZ = 'pixel index Z'
        else:
            valuesZ = None
            labelZ = None
        #We do not need to handle the qslice_viewer_option separately here since there is only one option right now
        if qspectrum_viewer_option == 0 and qslice_viewer_option == 0:  #Loadings
            mzSpectra = analysis_object['npg_labels_medianmz'][:]
            labelSpectra = "m/z"
            mzSlice = None
            labelSlice = None
        elif qspectrum_viewer_option > 0 and qslice_viewer_option > 0:
            mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ =\
                super(omsi_npg, cls).v_qmz(analysis_object, qslice_viewer_option - 1, qspectrum_viewer_option - 1)
        elif qspectrum_viewer_option == 0 and qslice_viewer_option > 0:
            mzSpectra = analysis_object['npg_labels_medianmz'][:]
            labelSpectra = "m/z"
            tempA, tempB, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ = \
                super(omsi_npg, cls).v_qmz(analysis_object, 0, qspectrum_viewer_option - 1)
        elif qspectrum_viewer_option > 0 and qslice_viewer_option == 0:
            mzSlice = analysis_object['npg_labels_medianmz'][:]
            labelSlice = "m/z"
            # Ignore the spatial axes. We need to use the axes of omsi_npg as defined above
            mzSpectra, labelSpectra, tempA, tempB, vX, lX, vY, lY, vZ, lZ = \
                super(omsi_npg, cls).v_qmz(analysis_object, 0, qspectrum_viewer_option - 1)

        return mzSpectra, labelSpectra, mzSlice, labelSlice, valuesX, labelX, valuesY, labelY, valuesZ, labelZ

    @classmethod
    def v_qspectrum_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_npg, cls).v_qspectrum_viewer_options(analysis_object)
        re = ["Peak array"] + dependent_options
        return re

    @classmethod
    def v_qslice_viewer_options(cls, analysis_object):
        """Define which viewer_options are supported for qspectrum URL's"""
        dependent_options = super(omsi_npg, cls).v_qslice_viewer_options(analysis_object)
        re = ["Peak array"] + dependent_options
        return re

    # ------------------------ viewer functions end ------------------------

    # main NPG function
    def execute_analysis(self):

         #Set default parameter values if needed
        if not self['npg_cluster_treecut']:
            self['npg_cluster_treecut'] = 0.1
        if not self['npg_mz_threshold']:
            self['npg_mz_threshold'] = 0.05
        if not self['npg_fasterhc_flag']:
            self['npg_fasterhc_flag'] = 3
        if not self['npg_cluster_method']:
            self['npg_cluster_method'] = 'median'
        if not self['npg_split_max'] and (fasterHC == 1):
            self['npg_split_max'] = 1000

        #Copy parameters to local variables for convenience
        peaksBins = self['peaksBins']
        peaksIntensities = self['npg_peaks_Intensities']
        peaksArrayIndex = self['npg_peaks_ArrayIndex']
        peaksMZData = self['peaksMZdata']
        peaksMZ = self['peaksMZ']
        MZ_TH = self['npg_mz_threshold']
        clusterCut = self['npg_cluster_treecut']
        fasterHC = self['npg_fasterhc_flag'] # set to 1 to perform split fastcluster, 2 for normal fastcluster, 3 for myHC
        SplitMax = self['npg_split_max']  # size limiter for HC clutering (smaller is usually faster) used if fasterHC == 1
        clusterMethod = unicode(self['npg_cluster_method'])  # Hierarchical Clustering agglomeration method

        # get dataset dimensions info from last peaksArrayIndex
        pAILast = peaksArrayIndex.shape[0] - 1
        Nx = peaksArrayIndex[pAILast][0] + 1
        Ny = peaksArrayIndex[pAILast][1] + 1

        print "[NPG]"
        print "NPG Parameters:"
        print "MZ_TH =", MZ_TH, "clusterCut =", clusterCut
        sys.stdout.flush()

        print "Performing first pass..."
        # each peak has a label to represent its cluster
        # all peaks are initialized to 0 = no label
        # regioncounter-1 is the last assigned cluster label
        peaksLabels = np.zeros_like(peaksMZ)
        peaksLabelsIndex = 0
        regioncounter = 1
        totalcounter = 0
        equivalentLabels = []
        LabelsList = []

        # pixelMap assigns a flag to each pixel where:
        # 0 = initialized value
        # 1 = peaks found by LPF
        # 2 = no peaks found by LPF
        # this allows us to skip pixels with no peaks both
        # when checking a pixel itself and its neighbors
        pixelMap = self.getPixelMap(Nx, Ny)

        timekeep1 = time()
        percentcheck = None
        # scan top-to-bottom
        for y in xrange(Ny - 1, -1, -1):
            # scan left-to-right
            for x in xrange(Nx):
                # print "pos: ", x, y
                percent = int(100. * float(totalcounter + 1) / float(Nx * Ny))
                timer = str(int(time() - timekeep1))
                eqs = str(len(equivalentLabels))
                regions = str(regioncounter)
                if (percent != percentcheck):
                    print "[", percent, "% -", timer, "s - eqs:", eqs, "- regs:", regions, "] \r",
                    sys.stdout.flush()
                    percentcheck = percent
                totalcounter += 1

                # check pixel map for pixels without peaks
                # if it doesnt have peaks, continue to next iteration
                if (pixelMap[x, y] != 1):
                    continue

                # get peaks for current coordinates
                myPeaks = self.getCoordPeaksB(x, y)

                # ----- check for 8-connected coord boundaries
                # check current pixel spatial position to account for boundaries
                CheckWest = CheckNorth = 1
                CheckNorthWest = CheckNorthEast = 1
                # dont check West coords on x = 0 boundary
                if (x == 0):
                    CheckWest = 0
                    CheckNorthWest = 0
                # dont check North coords on y = top boundary
                if (y == Ny - 1):
                    CheckNorth = 0
                    CheckNorthWest = 0
                    CheckNorthEast = 0
                # dont check East coords on x = top boundary
                if (x == Nx - 1):
                    CheckNorthEast = 0


                # prevent checking empty pixels
                if (CheckWest):
                    if (pixelMap[x - 1, y] != 1):
                        CheckWest = 0
                if (CheckNorth):
                    if (pixelMap[x, y + 1] != 1):
                        CheckNorth = 0
                if (CheckNorthWest):
                    if (pixelMap[x - 1, y + 1] != 1):
                        CheckNorthWest = 0
                if (CheckNorthEast):
                    if (pixelMap[x + 1, y + 1] != 1):
                        CheckNorthEast = 0

                # if all neighbors aren't available, add all peaks in pixel as distinct labels
                # Example: First pixel
                # initialize cluster labels to top-left coord
                # if(x == 0 and y == Ny-1):
                if (CheckWest == CheckNorth == CheckNorthWest == CheckNorthEast == 0):
                    # print "NON"
                    peaksLabelsIndex = self.getCoordIdxB(x, y)
                    for i in xrange(myPeaks.shape[0]):
                        peaksLabels[peaksLabelsIndex] = regioncounter
                        LabelsList.append(Node(regioncounter))
                        self.MakeSet(LabelsList[-1])
                        peaksLabelsIndex += 1
                        regioncounter += 1

                # check neighbors
                else:
                    # check Neighbors start ----------
                    peaksLabelsIndex = self.getCoordIdxB(x, y)

                    # gather corresponding neighbors's peak and label data
                    if (CheckWest):
                        [myWestPeaks, myWestLabels] = self.getCoordInfoB(x - 1, y, peaksLabels)
                    if (CheckNorth):
                        [myNorthPeaks, myNorthLabels] = self.getCoordInfoB(x, y + 1, peaksLabels)
                    if (CheckNorthWest):
                        [myNorthWestPeaks, myNorthWestLabels] = self.getCoordInfoB(x - 1, y + 1, peaksLabels)
                    if (CheckNorthEast):
                        [myNorthEastPeaks, myNorthEastLabels] = self.getCoordInfoB(x + 1, y + 1, peaksLabels)

                    # for each peak in current pixel, check the
                    # appropiate neighbor pixels for similar peaks
                    for i in xrange(myPeaks.shape[0]):
                        currentPeak = myPeaks[i]

                        westFlag = 0
                        northFlag = 0
                        northwestFlag = 0
                        northeastFlag = 0
                        newFlag = 0
                        nearLabels = []

                        # check west neighbor for similar peak
                        if (CheckWest):

                            # NOTE: getNearestPeakIndex finds the first index that is nearest
                            WIdx = self.getNearestPeakIndex(myWestPeaks, currentPeak)
                            Wpeak = myWestPeaks[WIdx]
                            WpeakLabel = myWestLabels[WIdx]

                            if (Wpeak >= currentPeak):
                                WThr = Wpeak - MZ_TH
                                if (WThr <= currentPeak):
                                    westFlag = 1
                                    nearLabels.append(WpeakLabel)
                            else:
                                WThr = Wpeak + MZ_TH
                                if (WThr >= currentPeak):
                                    westFlag = 1
                                    nearLabels.append(WpeakLabel)

                        # check north neighbor for similar peak
                        if (CheckNorth):
                            NIdx = self.getNearestPeakIndex(myNorthPeaks, currentPeak)
                            Npeak = myNorthPeaks[NIdx]
                            NpeakLabel = myNorthLabels[NIdx]

                            if (Npeak >= currentPeak):
                                NThr = Npeak - MZ_TH
                                if (NThr <= currentPeak):
                                    northFlag = 1
                                    nearLabels.append(NpeakLabel)
                            else:
                                NThr = Npeak + MZ_TH
                                if (NThr >= currentPeak):
                                    northFlag = 1
                                    nearLabels.append(NpeakLabel)

                        # check northwest neighbor for similar peak
                        if (CheckNorthWest):
                            NWIdx = self.getNearestPeakIndex(myNorthWestPeaks, currentPeak)
                            NWpeak = myNorthWestPeaks[NWIdx]
                            NWpeakLabel = myNorthWestLabels[NWIdx]

                            if (NWpeak >= currentPeak):
                                NWThr = NWpeak - MZ_TH
                                if (NWThr <= currentPeak):
                                    northwestFlag = 1
                                    nearLabels.append(NWpeakLabel)
                            else:
                                NWThr = NWpeak + MZ_TH
                                if (NWThr >= currentPeak):
                                    northwestFlag = 1
                                    nearLabels.append(NWpeakLabel)

                        # check northeast neighbor for similar peak
                        if (CheckNorthEast):
                            NEIdx = self.getNearestPeakIndex(myNorthEastPeaks, currentPeak)
                            NEpeak = myNorthEastPeaks[NEIdx]
                            NEpeakLabel = myNorthEastLabels[NEIdx]

                            if (NEpeak >= currentPeak):
                                NEThr = NEpeak - MZ_TH
                                if (NEThr <= currentPeak):
                                    northeastFlag = 1
                                    nearLabels.append(NEpeakLabel)
                            else:
                                NEThr = NEpeak + MZ_TH
                                if (NEThr >= currentPeak):
                                    northeastFlag = 1
                                    nearLabels.append(NEpeakLabel)


                        # ------- decide current peak label according to neighbor info
                        nearLabelsLen = len(nearLabels)
                        if (nearLabelsLen == 0):
                            # unique peak, add new label
                            choosenLabel = regioncounter
                            LabelsList.append(Node(regioncounter))
                            self.MakeSet(LabelsList[-1])
                            regioncounter += 1
                            newFlag = 1

                        elif (nearLabelsLen == 1):
                            # similar peak in only one 8-connected coord, assign it as label
                            choosenLabel = nearLabels[0]
                        elif (nearLabelsLen > 1):
                            # similar peak in several coordinates
                            # assign lower label to current peak
                            choosenLabel = min(nearLabels)

                            # all similar coordinates belong to same region
                            for t in xrange(nearLabelsLen - 1):
                                labelA = nearLabels[t]
                                labelB = nearLabels[t + 1]
                                if (labelA != labelB):
                                    # if two labels are different, union them to
                                    # note that both labels share same region
                                    # add relationship tupple to equivalentLabels
                                    equivalentLabels.append([labelA, labelB])
                                    self.Union(LabelsList[int(labelA - 1)], LabelsList[int(labelB - 1)])

                        else:
                            # this condition should never happen
                            # unless nearLabelsLen is negative
                            print "Error finding nearby labels"

                        # current peak label resolved, move index to next peak
                        peaksLabels[peaksLabelsIndex] = choosenLabel
                        peaksLabelsIndex += 1

                    # check Neighbors end ----------
                    # end of  [for i in xrange(myPeaks.shape[0]) ] ----------
                    # end of [else (if not init coord)] -------
                    # end of [for x in xrange(Nx)] -------
        # end of [for y in xrange(Ny-1,-1,-1)] -------
        print "\ndone!"
        sys.stdout.flush()

        # --- Starting Second Pass ---
        print "Performing second pass..."
        totalcounter = 0
        percentcheck = None
        timekeep1 = time()
        for i in xrange(peaksLabels.shape[0]):
            percent = int(100. * float(totalcounter + 1) / float(peaksLabels.shape[0]))
            timer = str(int(time() - timekeep1))
            if (percent != percentcheck):
                print "[", percent, "% -", timer, "s ]\r",
                sys.stdout.flush()
                percentcheck = percent

            # second pass replaces each label with root labels of union-find
            # this fixes the equivalences noted in the first pass
            mylabel = peaksLabels[i]
            mylabelIdx = int(mylabel - 1)
            peaksLabels[i] = int(str(self.Find(LabelsList[mylabelIdx])))

            totalcounter += 1
        print "\ndone!"
        sys.stdout.flush()

        # total regions/labels after resolving equivalences
        UniqueLabels = [int(str(self.Find(i))) for i in LabelsList]
        UniqueLabels = np.unique(UniqueLabels)
        print "# of regions:", UniqueLabels.shape[0]
        sys.stdout.flush()

        # -------- Start Hierarchical Clustering --------
        # -------- gather median info --------
        print "\n[HC]"
        print "Gathering NPG clusters info..."

        LabelsMedianMZ = self.getClustersInfo(peaksLabels, UniqueLabels)

        print "\ndone!"

        # ------------------------------------
        # --- Hierarchical Clustering (HC) ---
        print "Region Labels:", len(LabelsMedianMZ)
        print "Performing Hierarchical Clustering..."
        sys.stdout.flush()
        timekeeper = time()

        # Use custom hierarchical clustering implementation
        if (fasterHC == 3):
            TreeCut = self.myHC(LabelsMedianMZ, clusterCut)
            TreeCut = TreeCut - 1  # start labels from 0 instead of 1 temporarily

        # Use fastcluster module for hierarchical clustering
        elif (fasterHC == 2):
            LMData = LabelsMedianMZ
            LMData = np.reshape(LMData, (LMData.shape[0], 1))
            LinkageMatrix = fastcluster.linkage_vector(X=LMData, method=clusterMethod, metric='euclidean')
            TreeCut = hier.fcluster(LinkageMatrix, t=clusterCut, criterion='distance') - 1

        # Use fastcluster module for hierarchical clustering
        # but split the work into smaller parts using the clusterCut
        else:
            # split HC
            LMArgS = np.argsort(LabelsMedianMZ)
            LMData = LabelsMedianMZ[LMArgS]
            LMDataList = self.splitLabelsList(LMData, clusterCut, SplitMax)

            print "Splits:", len(LMDataList)
            sys.stdout.flush()

            PrevTCMax = 0
            FullTC = []
            percentcheck = None
            totalcounter = 0
            timekeep1 = time()
            for myLMelem in LMDataList:
                percent = int(100. * float(totalcounter + 1) / float(len(LMDataList)))
                timer = str(int(time() - timekeep1))
                if (percent != percentcheck):
                    print "[", percent, "% -", timer, "s ]\r",
                    sys.stdout.flush()
                    percentcheck = percent

                LMData = np.reshape(myLMelem, (myLMelem.shape[0], 1))
                LinkageMatrix = fastcluster.linkage_vector(X=LMData, method=clusterMethod, metric='euclidean')
                TreeCut = hier.fcluster(LinkageMatrix, t=clusterCut, criterion='distance')
                TreeCut += PrevTCMax
                PrevTCMax = TreeCut.max()
                if len(FullTC) == 0:
                    FullTC = TreeCut
                else:
                    FullTC = np.append(FullTC, TreeCut)

                totalcounter += 1

            LMArgS2 = np.argsort(LMArgS)  # reverse the sort
            TreeCut = FullTC[LMArgS2] - 1  # start labels from 0 instead of 1 temporarily

        print "\ndone! [", time() - timekeeper, "s ]"
        print "Reassigning labels..."
        sys.stdout.flush()

        peaksLabelsInt = peaksLabels.astype('int32')
        LabelsRange = np.arange(0, peaksLabelsInt.max() + 1)
        LabelsRange[UniqueLabels] = TreeCut
        HCpeaksLabels = LabelsRange[peaksLabelsInt] + 1  # make labels from 1 onward, leave 0 for space filling
        HCLabelsList = np.arange(1, HCpeaksLabels.max() + 1)
        HCpeaksLabels = HCpeaksLabels.astype('float32')

        print "done!"

        print "Collecting data into HDF5..."

        #Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
        #We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        #handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        #to ensure a consitent behavior we convert the values directly here


        # Results  -------
        # NPG Results
        self['npg_labels_medianmz'] = np.asarray(LabelsMedianMZ)
        self['npg_labels_list'] = np.asarray(UniqueLabels)
        self['npg_peaks_labels'] = np.asarray(peaksLabels)

        # NPG-HC Results
        # if(fasterHC == 2):
        # 	npghcLM = np.asarray( LinkageMatrix )
        # 	self.add_analysis_data( name='npghc_linkage_matrix' , data=npghcLM , dtype=str(npghcLM.dtype) )
        self['npghc_tc_labels'] = np.asarray(TreeCut + 1)
        self['npghc_labels_list'] = np.asarray(HCLabelsList)
        self['npghc_peaks_labels'] = np.asarray(HCpeaksLabels)

        print "Collecting done."
        print "--- finished ---"

    # get image slice z from peak data arrays
    def record_execute_analysis_outputs(self, analysis_output):
        """We are recording our outputs manually as part of the execute function"""
        pass


    # +++++++++++++++++++++++ helper functions ++++++++++++++++++++++++++
    @classmethod
    # get image slice z from peak data arrays
    def getnpgimage(cls, PeaksLabels, LabelsList, peaksArrayIndex, peaksIntensities, z):

        LabelID = LabelsList[z]

        # overall data
        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]
        Nx = peaksArrayIndex[pAILast][0] + 1
        Ny = peaksArrayIndex[pAILast][1] + 1

        img = np.zeros([Nx, Ny])

        mycounter = 0
        for xCoord in xrange(Nx):
            for yCoord in xrange(Ny):

                if (xCoord == pAImaxX and yCoord == pAImaxY):
                    pStart = peaksArrayIndex[pAILast][2]
                    pEnd = peaksIntensities.shape[0]
                else:
                    pAIXLeftIdx = xCoord * Ny
                    pAIindex = pAIXLeftIdx + yCoord
                    pStart = peaksArrayIndex[pAIindex][2]
                    pEnd = peaksArrayIndex[pAIindex + 1][2]

                rInts = peaksIntensities[pStart:pEnd]
                rLabels = PeaksLabels[pStart:pEnd]
                imgval = 0.0

                lblidx = np.where(rLabels == LabelID)[0]

                if (len(lblidx) == 0):
                    imgval = 0.0
                elif (len(lblidx) == 1):
                    imgval = rInts[lblidx[0]]
                else:
                    valInts = rInts[lblidx]
                    imgval = valInts.max()

                img[xCoord, yCoord] = imgval

        return img

    # get spectrum of an x/y coordinate from peak data arrays
    @classmethod
    def getnpgspec(cls, PeaksLabels, LabelsList, peaksArrayIndex, peaksIntensities, xCoord, yCoord):

        spec = np.zeros([len(LabelsList)])

        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]

        pStart = 0
        pEnd = 1

        # if its the last coord
        if (xCoord == pAImaxX and yCoord == pAImaxY):
            pStart = peaksArrayIndex[pAILast][2]
            pEnd = peaksIntensities.shape[0]

        # if its not last coord
        else:
            pAIXLeftIdx = np.searchsorted(peaksArrayIndex[:, 0], xCoord)
            pAIindex = pAIXLeftIdx + yCoord
            pStart = peaksArrayIndex[pAIindex][2]
            pEnd = peaksArrayIndex[pAIindex + 1][2]

        rLabels = PeaksLabels[pStart:pEnd]
        rInts = peaksIntensities[pStart:pEnd]

        for i in xrange(len(rLabels)):
            mylbl = rLabels[i]
            spec[LabelsList == mylbl] = rInts[i]

        return spec

    # custom hierarchical clustering
    def myHC(self, labelsMMz, TreeCut):
        thelist, inverse_idxs = np.unique(labelsMMz, return_inverse=True)
        lbls = np.arange(1, len(thelist) + 1)

        ltmz = np.copy(thelist)
        ltid = np.copy(lbls)

        dif = np.diff(thelist)

        print "NPG Amount : ", len(thelist)

        percentcheck = None
        totalcounter = 0
        totallen = len(thelist)
        timekeep1 = time()
        while (len(thelist) > 1):

            # -- counter + --
            percent = int(100. * float(totalcounter) / float(totallen))
            timer = int(time() - timekeep1)
            if (percent != percentcheck):
                print "[", percent, "% - dif.min:", round(dif.min(), 4), "len:", len(thelist), "\t-", timer, "s ]\r",
                sys.stdout.flush()
                percentcheck = percent
            totalcounter += 1
            # -- counter - --

            # stop if treecut is reached
            if (dif.min() >= TreeCut):
                break

            mpos = dif.argmin()
            thelistL = thelist[mpos]
            thelistR = thelist[mpos + 1]

            # replace the 2 nearest with median of them
            mymedian = np.median([thelistL, thelistR])
            thelist[mpos] = mymedian
            thelist = np.delete(thelist, mpos + 1)

            # update dif
            dif = np.delete(dif, mpos)
            # update dif left
            if (mpos != 0):
                dif[mpos - 1] = thelist[mpos] - thelist[mpos - 1]
            # update dif right
            if (mpos != len(dif)):
                dif[mpos] = thelist[mpos + 1] - thelist[mpos]

            lblsL = lbls[mpos]
            lblsR = lbls[mpos + 1]

            mmin = min(lbls[mpos], lbls[mpos + 1])
            mmax = max(lbls[mpos], lbls[mpos + 1])

            # re-label left position as lower label
            lbls[mpos] = mmin
            # mark previous higher labels as lower labels
            idxs = np.where(lbls == mmax)[0]
            lbls[idxs] = mmin
            # remove the higher label
            lbls = np.delete(lbls, mpos + 1)

            idxs = np.where(ltid == mmax)[0]
            ltid[idxs] = mmin

        # ---- end of while

        # --- labeling ---
        peaksLabels = ltid[inverse_idxs]
        UniqueLabels = np.unique(peaksLabels)

        ClusterRange = np.arange(1, len(UniqueLabels) + 1)
        LabelsRange = np.arange(0, peaksLabels.max() + 1)
        LabelsRange[UniqueLabels] = ClusterRange
        result = LabelsRange[peaksLabels]

        return result

    # split Labels into ammounts limited by SplitMax
    def splitLabelsList(self, LabelData, Thres, SplitMax):

        myLabels = np.copy(LabelData)
        myLabelsList = []
        stopFlag = 0
        snum = 0

        try:
            snum = np.ceil(len(myLabels) / float(SplitMax))
            myAmnt = int(len(myLabels) / snum)
        except:
            print "\nError splitting labels list with SplitMax =", SplitMax
            stopFlag = 1
            myLabelsList.append(myLabels)

        while ( (len(myLabelsList) < snum) and (stopFlag == 0) ):
            try:
                RightSide = myLabels[myAmnt:]
                myDiff = np.diff(RightSide)
                SplitPos = myAmnt + (np.where(myDiff > Thres)[0][0] + 1)
                myLabelsList.append(myLabels[0:SplitPos])
                myLabels = myLabels[SplitPos:]
            except:
                # if something goes wrong, dont split array and stop further splits
                stopFlag = 1
                myLabelsList.append(myLabels)

        return myLabelsList

    # gathers information of clusters and saves the data
    def getClustersInfo(self, GpeaksLabels, GLabelsList):
        import bisect

        peaksBins = self['peaksBins']
        peaksMZdata = self['peaksMZdata']

        labelsArgs = np.argsort(GpeaksLabels)
        SortedPL = GpeaksLabels[labelsArgs]
        LabelsMedianMZ = np.empty_like(GLabelsList, dtype=float)
        idxLeft = 0

        timekeep = time()
        percentcheck = None
        for i in xrange(len(GLabelsList)):
            currentlabel = GLabelsList[i]

            idxRight = bisect.bisect_right(SortedPL, currentlabel)
            currentMZs = peaksMZdata[peaksBins[labelsArgs[idxLeft:idxRight]]]
            LabelsMedianMZ[i] = np.median(currentMZs)
            idxLeft = idxRight

            percent = int(100. * float(i + 1) / float(len(GLabelsList)))
            timer = int(time() - timekeep)
            if (percent != percentcheck):
                print "[", i + 1, "of", len(GLabelsList), "-", percent, "% -", timer, "s ] \r",
                sys.stdout.flush()
                percentcheck = percent

        return LabelsMedianMZ

    # find nearest element
    def getNearestPeakIndex(self, myPeaksArray, myPeak):
        nearestPeakIndex = (np.abs(myPeaksArray - myPeak)).argmin()
        return nearestPeakIndex

    # binary search
    # returns the peaks of a coordinate (x,y) from peaksArrayIndex
    def getCoordPeaksB(self, xCoord, yCoord):
        peaksArrayIndex = self['npg_peaks_ArrayIndex']
        peaksMZ = self['peaksMZ']

        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]

        pStart = 0
        pEnd = 1

        # check its a valid coordinate
        if (xCoord < 0 or xCoord > pAImaxX or yCoord < 0 or yCoord > pAImaxY):
            print "Error: Invalid Coordinate"
            return 0

        else:
            # if its the last coord
            if (xCoord == pAImaxX and yCoord == pAImaxY):
                pStart = peaksArrayIndex[pAILast][2]
                pEnd = peaksMZ.shape[0]

            # if its not last coord
            else:
                pAIXLeftIdx = np.searchsorted(peaksArrayIndex[:, 0], xCoord)
                pAIindex = pAIXLeftIdx + yCoord
                pStart = peaksArrayIndex[pAIindex][2]
                pEnd = peaksArrayIndex[pAIindex + 1][2]

            return peaksMZ[pStart:pEnd]

    # binary search
    # returns the peak index of a coordinate (x,y) from peaksArrayIndex
    def getCoordIdxB(self, xCoord, yCoord):
        peaksArrayIndex = self['npg_peaks_ArrayIndex']
        peaksMZ = self['peaksMZ']

        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]

        pStart = 0
        pEnd = 1

        # check its a valid coordinate
        if (xCoord < 0 or xCoord > pAImaxX or yCoord < 0 or yCoord > pAImaxY):
            print "Error: Invalid Coordinate"
            return 0

        else:
            # if its the last coord
            if (xCoord == pAImaxX and yCoord == pAImaxY):
                pStart = peaksArrayIndex[pAILast][2]
                pEnd = peaksMZ.shape[0]

            # if its not last coord
            else:
                pAIXLeftIdx = np.searchsorted(peaksArrayIndex[:, 0], xCoord)
                pAIindex = pAIXLeftIdx + yCoord
                pStart = peaksArrayIndex[pAIindex][2]

            return pStart

    # binary search
    # 1. rPeaks = returns the peaks of a coordinate (x,y) from peaksArrayIndex
    # 2. rLabels = returns the labels of the rPeaks
    def getCoordInfoB(self, xCoord, yCoord, peaksLabels):
        peaksArrayIndex = self['npg_peaks_ArrayIndex']
        peaksMZ = self['peaksMZ']

        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]

        pStart = 0
        pEnd = 1

        # check its a valid coordinate
        if (xCoord < 0 or xCoord > pAImaxX or yCoord < 0 or yCoord > pAImaxY):
            print "Error: Invalid Coordinate"
            return 0

        else:
            # if its the last coord
            if (xCoord == pAImaxX and yCoord == pAImaxY):
                pStart = peaksArrayIndex[pAILast][2]
                pEnd = peaksMZ.shape[0]

            # if its not last coord
            else:
                pAIXLeftIdx = np.searchsorted(peaksArrayIndex[:, 0], xCoord)
                pAIindex = pAIXLeftIdx + yCoord
                pStart = peaksArrayIndex[pAIindex][2]
                pEnd = peaksArrayIndex[pAIindex + 1][2]

            rPeaks = peaksMZ[pStart:pEnd]
            rLabels = peaksLabels[pStart:pEnd]

            return rPeaks, rLabels

    # pixelMap assigns a flag to each pixel where:
    # 0 = initialized value
    # 1 = peaks found by LPF
    # 2 = no peaks found by LPF
    # 3 = error
    # this allows us to skip pixels with no peaks both
    # when checking a pixel itself and its neighbors
    def getPixelMap(self, Nx, Ny):
        peaksArrayIndex = self['npg_peaks_ArrayIndex']
        peaksMZ = self['peaksMZ']

        pixelMap = np.zeros((Nx, Ny))

        pAILast = peaksArrayIndex.shape[0] - 1
        pAImaxX = peaksArrayIndex[pAILast][0]
        pAImaxY = peaksArrayIndex[pAILast][1]

        pStart = 0
        pEnd = 0
        pAIindex = 0

        for x in xrange(Nx):
            for y in xrange(Ny):
                # pStart is the index location of
                # current pixel's peaks values
                pStart = peaksArrayIndex[pAIindex][2]

                # pEnd is the index location of next pixel's values
                # if current pixel its the last element of peaksArrayIndex
                # we set pEnd peaksMZ size to prevent out of bounds
                if (x == pAImaxX and y == pAImaxY):
                    pEnd = peaksMZ.shape[0]
                else:
                    pEnd = peaksArrayIndex[pAIindex + 1][2]

                # Compare pStart and pEnd to obtain the
                # ammount of peaks in current pixel and
                # set the appropiate value to pixelMap
                if (pEnd > pStart):
                    # some peaks in current pixel
                    pixelMap[x, y] = 1
                elif (pEnd == pStart):
                    # no peaks in current pixel
                    pixelMap[x, y] = 2
                else:
                    # error
                    pixelMap[x, y] = 3

                pAIindex += 1

        return pixelMap

    # union find
    # http://code.activestate.com/recipes/577225-union-find/

    def MakeSet(self, x):
        x.parent = x
        x.rank = 0

    def Union(self, x, y):
        xRoot = self.Find(x)
        yRoot = self.Find(y)
        if xRoot.rank > yRoot.rank:
            yRoot.parent = xRoot
        elif xRoot.rank < yRoot.rank:
            xRoot.parent = yRoot
        elif xRoot != yRoot:  # Unless x and y are already in same set, merge them
            yRoot.parent = xRoot
            xRoot.rank = xRoot.rank + 1

    def Find(self, x):
        if x.parent == x:
            return x
        else:
            x.parent = self.Find(x.parent)
            return x.parent


def main(argv=None):
    if argv is None:
        argv = sys.argv

        # Check for correct usage
    if len(argv) != 5:
        print "USAGE: Call \"omsi_npg OMSI_FILE [expIndex dataIndex LPFanalysisIndex]   \" "
        print "Nearby-Peaks Global peakfinding"
        exit(0)

    #Read the input arguments
    omsiInFile = argv[1]
    expIndex = int(argv[2])
    dataIndex = int(argv[3])
    analysisIndex = int(argv[4])

    #Open the input HDF5 file
    try:
        omsiFile = omsi_file(omsiInFile)
    except:
        print "Error opening input file:", sys.exc_info()[0]
        exit(0)

    print "Input file: ", omsiInFile
    print "LPF Analysis Index: ", analysisIndex

    #Get the experiment and data
    exp = omsiFile.get_experiment(expIndex)
    data = exp.get_msidata(dataIndex)
    analysis = exp.get_analysis(analysisIndex)

    print "\n[loading data...]"
    peaksBins = analysis['LPF_Peaks_MZ'][:]
    peaksIntensities = analysis['LPF_Peaks_Vals'][:]
    peaksArrayIndex = analysis['LPF_Peaks_ArrayIndex'][:]
    peaksMZdata = data.mz[:]
    peaksMZ = peaksMZdata[peaksBins]

    print "[done!] lpf data shapes:"
    print "peaksBins shape: ", peaksBins.shape
    print "peaksIntensities shape: ", peaksIntensities.shape
    print "peaksArrayIndex shape: ", peaksArrayIndex.shape
    print "peaksMZdata shape: ", peaksMZdata.shape
    print "peaksMZ shape: ", peaksMZ.shape
    sys.stdout.flush()

    # npg
    myNPG = omsi_npg(name_key="omsi_npg_" + str(ctime()))
    print "--- Executing NPG ---"
    #myNPG.omsi_npg_exec(peaksBins, peaksIntensities, peaksArrayIndex, peaksMZdata, peaksMZ)
    myNPG.execute(peaksBins=peaksBins,
                  npg_peaks_Intensities=peaksIntensities,
                  npg_peaks_ArrayIndex=peaksArrayIndex,
                  peaksMZdata=peaksMZdata,
                  peaksMZ=peaksMZ)
    print "\nResults:"
    NPGPL = myNPG['npghc_peaks_labels']
    print "NPG HC Peaks Labels: \n", NPGPL
    NPGLL = myNPG['npghc_labels_list']
    print "NPG HC Labels List: \n", NPGLL

    print "\nsaving HDF5 analysis..."
    NPGanalysis, analysisindex = exp.create_analysis(myNPG)
    print "done!"

    print "npg analysis index:", analysisindex
    print "omsi_npg complete for input file: ", omsiInFile, "\n"


# stop python, used for debugging
def stop():
    raw_input("Stop!")


if __name__ == "__main__":
    main()
