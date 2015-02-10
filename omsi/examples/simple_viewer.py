"""Simple viewer for OpenMSI data"""

from omsi.dataformat.omsi_file.main_file import omsi_file
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


class MyViewer : 
    def __init__(self,  data , mzdata ) :
        """Create a simple viewer with image of the data and a cruve plot for a spectrum.
           
           :param data: Reference to the hdf5 dataset of the image
           :param mzdata: mz values of the instrument to be displayed as axis in the curve plot. \
                          may be None in case the mz data is unknown
        """
        #Compute the reference image using the full dataset
        print "Computing max projection of all spectra"
        maxVal = np.max( data[:,:,:] , axis=2 )
        print "Initalizing plots"
        #Create plot of the physical image
        self.hdfdata = data
        self.mzdata = mzdata
        self.mainFig = plt.figure()
        gs = gridspec.GridSpec(1, 2,width_ratios=[1,2])
        self.imageFig = self.mainFig.add_subplot(gs[0])
        self.imageFig.autoscale(True,'both',tight=True)
        self.imagePlot = self.imageFig.pcolor( np.log( maxVal ) ) 
        #Create curve plot for a single spectrum         
        self.curveFig = self.mainFig.add_subplot(gs[1])
        self.curve = self.curveFig.plot( data[0,0,:]   )[0]            
        self.curve.set_ydata( data[0,0,:] )
        if mzdata is not None :
            self.curve.set_xdata( self.mzdata )
            self.curveFig.set_xlim( np.min( self.mzdata ) , np.max( self.mzdata) )
            self.curveFig.set_xlabel("m/z")
        self.imageFig.figure.canvas.mpl_connect( 'button_press_event' , self)
        print "Done creating the view"
   
    def __call__(self, event ):   
        """Callback function used to update the curve plot when clicking on the image plot"""
        xloc = int(event.ydata) #The x and y axis appear to be switched in the data and image
        yloc = int( event.xdata) 
        if xloc < self.hdfdata.shape[0] and yloc < self.hdfdata.shape[1] :
            #Load the data for the selected spectrum and update the curve plot
            newData = self.hdfdata[ int(event.ydata) ,int( event.xdata) , :  ]
            maxVal = np.max( newData )
            self.curve.set_ydata( newData ) 
            self.curveFig.set_ylim(  0 , maxVal )
            self.curve.figure.canvas.draw()

def main(argv=None):
    """Then main function"""

    import sys
    from sys import argv,exit
    
    if argv is None:
        argv = sys.argv
   
    #Check for correct usage
    if len(argv) <2 :
        print "USAGE: Call \"omsiHDF5File [expIndex dataIndex]   \" "
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
    exp = omsiFile.get_experiment( expIndex )
    data = exp.get_msidata(dataIndex)
    mzdata = exp.get_instrument_info().get_instrument_mz()    

    if data is None:
        print "Could not access image data for the experiment"
        exit(0)

    #Initalize the viewer
    viewer = MyViewer(  data , mzdata )
    plt.show()     




if __name__ == "__main__":
    main() 

