from nmf import *
from numpy import *
from omsi.analysis.omsi_analysis_base import omsi_analysis_base
from omsi.analysis.omsi_analysis_data import omsi_analysis_data
from omsi.shared.omsi_dependency import *
from omsi.dataformat.omsi_file import *

class omsi_nmf(omsi_analysis_base) :
    """Class defining a basic nmf analysis for a 2D MSI data file or slice of the data"""

    def __init__(self, nameKey="undefined"):
        """Initalize the basic data members"""
        super(omsi_nmf,self).__init__()
        self.analysis_identifier = nameKey
        
    @classmethod
    def v_qslice(cls , anaObj , z , viewerOption=0) :
        """Implement support for qslice URL requests for the viewer""" 
        if viewerOption == 0 :
            dataset = anaObj[ 'ho' ]
            try :
                data = eval("dataset[:,:, %s]" %(z,))
                return data
            except:
                return None
        elif viewerOption > 0 :
            return super(omsi_nmf,cls).v_qslice( anaObj , z, viewerOption-1)
    
    @classmethod
    def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
        """Implement support for qspectrum URL requests for the viewer"""
        from omsi.shared.omsi_data_selection import check_selection_string, selection_type, selection_to_indexlist
        data = None
        customMZ = None
        if viewerOption == 0 : #Loadings
            dataset = anaObj[ 'ho' ]
            if (check_selection_string(x) == selection_type['indexlist']) and (check_selection_string(y) == selection_type['indexlist']) :
                data = eval("dataset[:][%s,%s, :]" %(x,y)) #Load the full data if multiple lists are given for the selection and let numpy handle the subselection
            else :
                data = eval("dataset[%s,%s, :]" %(x,y))
            customMZ = None
        elif viewerOption > 0 :
            return super(omsi_nmf,cls).v_qspectrum( anaObj , x , y, viewerOption-1)
            
        return data, customMZ
        
    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
        """Implement support for qmz URL requests for the viewer"""
        mzSpectra =  None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        #We do not need to handle the qslice_viewerOption separately here since there is only one option right now
        if qspectrum_viewerOption == 0 and qslice_viewerOption==0: #Loadings
            mzSpectra = arange( 0 ,  anaObj[ 'ho' ].shape[2] )
            labelSpectra = "Component Index"
            mzSlice  = None
            labelSlice = None
        elif qspectrum_viewerOption > 0 and qslice_viewerOption>0 :
            mzSpectra, labelSpectra, mzSlice, labelSlice = super(omsi_nmf,cls).v_qmz( anaObj, qslice_viewerOption-1 , qspectrum_viewerOption-1)
        elif qspectrum_viewerOption == 0 and qslice_viewerOption>0 :
            mzSpectra = arange( 0 ,  anaObj[ 'ho' ].shape[2] )
            labelSpectra = "Component Index"
            tempA, tempB, mzSlice, labelSlice = super(omsi_nmf,cls).v_qmz( anaObj, 0 , qspectrum_viewerOption-1)
        elif qspectrum_viewerOption > 0 and qslice_viewerOption==0 :
            mzSlice =  arange( 0 ,  anaObj[ 'ho' ].shape[2] )
            labelSlice = "Component Index"
            mzSpectra, labelSpectra, tempA, tempB = super(omsi_nmf,cls).v_qmz( anaObj, 0 , qspectrum_viewerOption-1)
        
        return mzSpectra, labelSpectra, mzSlice, labelSlice
    
    @classmethod
    def v_qspectrum_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        dependent_options = super(omsi_nmf,cls).v_qspectrum_viewerOptions(anaObj)
        re = ["NMF Loadings"] + dependent_options
        return re
        
    @classmethod
    def v_qslice_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        dependent_options = super(omsi_nmf,cls).v_qslice_viewerOptions(anaObj)
        re = ["NMF Images"] + dependent_options
        return re
        
        
    def execute_nmf(self, msidata , numComponents=20, timeOut=600, numIter=2000, tolerance=0.0001, msidata_dependency=None) :
        """Execute the nmf for the given msidata
        
           Keyword Arguments:
           :param msidata: h5py data object with the msi data.
           :param numComponents: Number of components to compute
           :param timeOut: Maximum time it should take
           :param numIter: Number of iterations
           :param tolerance: tolerance for a relative stopping condition
           :param msidata_dependency: The h5py.Datase or the omsi_file_analysis or omsi_file_msidata object on which the msidata input data object depends/originates from. 
           
        """
        #Copy the input data
        data = msidata[:]
        nx= data.shape[0]
        ny= data.shape[1]
        #Determine the input shape
        numBins   = data.shape[-1]
        numPixels = data.size / numBins #The last data dimension is assumed to contain the spectra
        
        #Reshape the data         
        data = data.reshape(numPixels, numBins)
        data = data.transpose()
        temp = data.max(axis = 1)
        idx = argwhere(temp>1000)
        idx = ravel(idx)
        data = data[idx,:];
        Nz = idx.shape[0];

        #Execute nmf
        (wo,ho) = nmf(data, random.randn(Nz,numComponents), random.randn(numComponents ,numPixels), tolerance, timeOut, numIter)

        #Reshape the ho matrix to be a 3D image cube
        ho = ho.transpose() 
        ho = ho.reshape( (nx,ny,numComponents))
        
        #Clear any previously stored analysis data
        self.clear_analysis_data()
        self.clear_parameter_data()
        self.clear_dependency_data()
        
        #Save the analysis data to the __data_list so that the data can be saved automatically by the omsi HDF5 file API
        self.add_analysis_data( name='wo' , data=wo , dtype=str(wo.dtype) )
        self.add_analysis_data( name='ho' , data=ho , dtype=str(ho.dtype) )
        
        #Save the analysis parameters to the __parameter_list so that the data can be saved automatically by the omsi HDF5 file API
        nc_param = asarray( [numComponents] )
        self.add_parameter_data( name='numComponents' , data=nc_param , dtype=str(nc_param.dtype) )
        t_param = asarray( [timeOut] )
        self.add_parameter_data( name='timeOut' , data=t_param , dtype=str(t_param.dtype) )
        ni_param = asarray( [numIter] )
        self.add_parameter_data( name='numIter' , data=ni_param , dtype=str(ni_param.dtype) )
        to_param = asarray( [tolerance] )
        self.add_parameter_data( name='tolerance' , data=to_param , dtype=str(to_param.dtype) )
        
        #Save the analysis dependencies to the __dependency_list so that the data can be saved automatically by the omsi HDF5 file API
        if msidata_dependency is not None :
            if isinstance( msidata_dependency , omsi_dependency ) :
                #Add the depency as given by the user
                self.add_dependency_data( msidata_dependency )
            else :
                #The user only gave us the object that we depend on so we need to construct the 
                self.add_dependency_data( omsi_dependency( param_name = 'msidata', link_name='msidata', omsi_object=msidata_dependency, selection=None ) )

           

def main(argv=None):
    '''Then main function'''

    import sys
    from sys import argv,exit
    from omsi.dataformat.omsi_file import omsi_file
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    if argv is None:
        argv = sys.argv

    #Check for correct usage
    if len(argv) <2 :
        print "USAGE: Call \"omsi_nmf [expIndex dataIndex]   \" "
        print "\n"
        print "This is a simple viewer for looking at OMSI HDF5 files."
        print "The viewer takes the index of the experiment and the"
        print "dataset to be used as optional input"
        exit(0)

    #Read the input arguments 
    omsiOutFile  = argv[1]
    expIndex = 0
    dataIndex = 0
    if len(argv)==4 :
        expIndex = int(argv[2] )
        dataIndex = int(argv[3] )

    #Open the input HDF5 file
    try:
        omsiFile = omsi_file( omsiOutFile , 'r' ) #Open file in read only mode
    except:
        print "Unexpected error creating the output file:", sys.exc_info()[0]
        exit(0)

    #Get the experiment and dataset
    exp = omsiFile.get_exp( expIndex )
    data = exp.get_msidata( dataIndex )

    #Execute the nmf
    testNMF = omsi_nmf()
    print "Executing nmf analysis"
    testNMF.execute_nmf( data )
    print "Getting nmf analysis results"
    wo = testNMF.get_analysis_data_by_name( 'wo' )['data']
    ho = testNMF.get_analysis_data_by_name( 'ho' )['data']
    print ho
    print "Plotting nmf analysis results"
    Nx = data.shape[0]
    Ny = data.shape[1]
    ho1 = ho[0,:];
    ho1 = ravel(ho1)
    ho1 = ho1.reshape(Nx,Ny)
    
    ho2 = ho[1,:];
    ho2 = ravel(ho2)
    ho2 = ho2.reshape(Nx,Ny)

    ho3 = ho[2,:];
    ho3 = ravel(ho3)
    ho3 = ho3.reshape(Nx,Ny)

    mainFig = plt.figure()
    gs = gridspec.GridSpec(1, 4)
    imageFig = mainFig.add_subplot(gs[0])
#    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( log (ho1 + 1) ) 

    imageFig = mainFig.add_subplot(gs[1])
#    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( log (ho2 + 1) )

    imageFig = mainFig.add_subplot(gs[2])
#    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( log (ho3 + 1) )


#    do the three color
    ho = ho.transpose()
    ho = ho.reshape(Nx,Ny,3)
    temp = log ( ho[:,:,0] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,0] = temp
    
    temp = log ( ho[:,:,1] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,1] = temp
    
    temp = log ( ho[:,:,2] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,2] = temp
    
    imageFig = mainFig.add_subplot(gs[3])
    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.imshow( ho )
   
    plt.show()


if __name__ == "__main__":
    main()

