from omsi.dataformat.omsi_file.main_file import omsi_file
import numpy as np
try :
    from PIL import Image
except :
    import Image
import datetime
import os, sys

def main(argv=None):
    """Then main function"""

    import sys
    from sys import argv,exit

    if argv is None:
        argv = sys.argv

    #Check for correct usage
    if len(argv) <3 :
        printHelp()
        exit(0)

    infile = argv[1]
    expIndex = int(argv[2])
    anaIndex = int(argv[3])
    outfileName = argv[4]

    #Open the file and get to the data
    #Note the analysis and experiment index may vary for different datasets
    f = omsi_file( infile , 'r' )
    e = f.get_experiment(expIndex)
    ana = e.get_analysis(anaIndex)
    #Get the peak finding data abd load it all
    pc = ana[ 'peak_cube' ]
    pm = ana[ 'peak_mz' ]

    #Load all the peak finding data
    peakCube = pc[:]
    peakMZ   = pm[:]
    xdim = peakCube.shape[0]
    ydim = peakCube.shape[1]
    numSlices = peakCube.shape[2]

    #Generate the images
    print "Generating images"
    names = []
    for i in range(0,numSlices) :
        sys.stdout.write("[" +str( int( 100.* float(i)/float(numSlices) )) +"%]"+ "\r")
        sys.stdout.flush()

        imData = peakCube[:,:,i]
        imData = imData.reshape((xdim,ydim))
        imData = imData / np.max(imData)
        a = Image.fromarray( imData.astype('float') * 255 )
        name = outfileName + str(peakMZ[i])
        name = name.replace( "." , "_" ) #Replace the "." in the string to make sure pdflatex does not complain about the filename
        name = name +".png"
        names.append( name )
        a.convert('L').save( name , 'PNG' )

    #Generate the latex file
    print "Generating LaTeX file"

    texFileName = outfileName+"summary.tex"
    latexFile = open( texFileName  , 'w' )
    latexFile.write("\documentclass[a4paper,10pt]{article}\n")
    latexFile.write("\usepackage[utf8x]{inputenc}\n")
    latexFile.write("\usepackage{graphicx}\n")
    latexFile.write("\\title{Peak Images}\n")
    latexFile.write("\\author{Oliver Ruebel}\n")
    latexFile.write("\date{"+str(datetime.date.today())+"}\n")
    latexFile.write("\\begin{document}\n")
    latexFile.write("\maketitle\n ")

    #Write the latex table with the image files
    numRows = 6
    numCols = 3
    totalNumRows = int( len(names) / float(numCols) +0.5 )
    for i in range( 0 , totalNumRows ):
        sys.stdout.write("[" +str( int( 100.* float(i)/float(totalNumRows) )) +"%]"+ "\r")
        sys.stdout.flush()

        if i % numRows == 0 :
            latexFile.write("\\begin{center}\n")
            latexFile.write("\\begin{tabular}{lll}\n")

        for ci in xrange(0,numCols) :
            index = i*numCols+ci
            if index < peakMZ.size :
                latexFile.write(str(peakMZ[index]))
            else :
                latexFile.write(" ")
            if ci < (numCols-1) :
                latexFile.write(" & ")
            else :
                latexFile.write( " \\\\ \n")

        for ci in xrange(0,numCols) :
            index = i*numCols+ci
            if index < len(names) :
                latexFile.write("\includegraphics[width=0.3\\textwidth]{"+names[index]+"}")
            else :
                latexFile.write(" ")
            if ci < (numCols-1) :
                latexFile.write(" & ")
            else :
                latexFile.write( " \\\\ \n")

        if i%numRows == (numRows-1) or i == (totalNumRows-1) :
            latexFile.write("\\end{tabular}\n")
            latexFile.write("\\end{center}\n")

    latexFile.write("\\end{document}\n")

def printHelp() :

    print "USAGE: Call \"python peakCubeOverview.py HDF5File expIndex anaIndex outFile\" "
    print ""
    print "HDF5File : The OpenMSI HDF5 file to be used"
    print "expIndex : The index of the experiment within the file"
    print "anaIndex : The index of the analysis with the global peak finding results"
    print "outFile : The basename (and path) for the output files"
    print ""
    print "Example usage:"
    print "python peakCubeOverview.py /work2/bowen/2012_0403_KBL_platename_100212.h5 0 1 /work2/bowen/test_images/slice"
    print ""
    print "Output : "
    print "1) One image for each m/z slice of the peak_cube data"
    print "2) Latex file with a summary of all images. The latex file can be compiled to pdf using pdflatex."

if __name__ == "__main__":
    main()
