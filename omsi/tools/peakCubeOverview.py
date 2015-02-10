"""
Simple helper tool used to generate a set of PNG images for a global peak
analysis (one per global peak) as well as a LaTeX document that summarizes
all the images in a single document.

NOTE: The module will try to build the LaTeX document using pdflatex, i.e,
      it is assumed that pdflatex is available.
"""

import numpy as np

from omsi.dataformat.omsi_file.main_file import omsi_file

try:
    from PIL import Image
except ImportError:
    import Image
import datetime
import os
import subprocess
import sys


def main(argv=None):
    """Then main function"""

    if argv is None:
        argv = sys.argv

    # Check for correct usage
    if len(argv) < 3:
        print_help()
        sys.exit(0)

    infile = argv[1]
    experiment_index = int(argv[2])
    analysis_index = int(argv[3])
    output_filename = argv[4]

    # Open the file and get to the data
    # Note the analysis and experiment index may vary for different datasets
    omsi_input_file = omsi_file(infile, 'r')
    input_experiment = omsi_input_file.get_experiment(experiment_index)
    input_analysis = input_experiment.get_analysis(analysis_index)
    # Get the peak finding data abd load it all
    input_dataset = input_analysis['peak_cube']
    input_mz_dataset = input_analysis['peak_mz']
    print input_analysis.items()

    # Load all the peak finding data
    peak_cube = input_dataset[:]
    peak_mz = input_mz_dataset[:]
    x_dim = peak_cube.shape[0]
    y_dim = peak_cube.shape[1]
    num_slices = peak_cube.shape[2]

    # Generate the images
    print "Generating images"
    names = []
    for row_index in range(0, num_slices):
        sys.stdout.write("[" + str(int(100. * float(row_index)/float(num_slices))) + "%]" + "\r")
        sys.stdout.flush()
        image_data = peak_cube[:, :, row_index]
        image_data = image_data.reshape((x_dim, y_dim))
        image_data = image_data / np.max(image_data)
        output_images = Image.fromarray(image_data.astype('float') * 255)
        name = output_filename + str(peak_mz[row_index])
        # Replace the "." in the string to make sure pdflatex does not complain about the filename
        name = name.replace(".", "_")
        name += ".png"
        names.append(name)
        output_images.convert('L').save(name, 'PNG')

    # Generate the latex file
    print "Generating LaTeX file"

    tex_filename = output_filename+"summary.tex"
    latex_file = open(tex_filename, 'w')
    latex_file.write("\documentclass[a4paper,10pt]{article}\n")
    latex_file.write("\usepackage[utf8x]{inputenc}\n")
    latex_file.write("\usepackage{graphicx}\n")
    latex_file.write("\\title{Peak Images}\n")
    latex_file.write("\\author{Oliver Ruebel}\n")
    latex_file.write("\date{"+str(datetime.date.today())+"}\n")
    latex_file.write("\\begin{document}\n")
    latex_file.write("\maketitle\n ")

    # Write the latex table with the image files
    number_of_columns = 3
    res = float(x_dim) / float(y_dim)
    number_of_rows = int(1.29 / (0.3*(float(x_dim)/float(y_dim))))
    if number_of_rows < 1:
        number_of_rows = 1
        print "WARNING: Images may not fit on a single page. Reduce the width of image plots"
    total_number_of_rows = int(len(names) / float(number_of_columns) + 0.5)
    for row_index in range(0, total_number_of_rows):
        sys.stdout.write("[" + str(int(100. * float(row_index)/float(total_number_of_rows))) + "%]" + "\r")
        sys.stdout.flush()
        if row_index % number_of_rows == 0:
            latex_file.write("\\begin{center}\n")
            latex_file.write("\\begin{tabular}{lll}\n")

        for column_index in xrange(0, number_of_columns):
            index = row_index*number_of_columns+column_index
            if index < peak_mz.size:
                latex_file.write(str(peak_mz[index]))
            else:
                latex_file.write(" ")
            if column_index < (number_of_columns-1):
                latex_file.write(" & ")
            else:
                latex_file.write(" \\\\ \n")

        for column_index in xrange(0, number_of_columns):
            index = row_index*number_of_columns+column_index
            if index < len(names):
                latex_file.write("\includegraphics[width=0.3\\textwidth]{"+os.path.basename(names[index])+"}")
            else:
                latex_file.write(" ")
            if column_index < (number_of_columns-1):
                latex_file.write(" & ")
            else:
                latex_file.write(" \\\\ \n")

        if row_index % number_of_rows == (number_of_rows-1) or row_index == (total_number_of_rows-1):
            latex_file.write("\\end{tabular}\n")
            latex_file.write("\\end{center}\n")

    latex_file.write("\\end{document}\n")
    latex_file.close()

    print "Building the PDF file"
    abs_path = os.path.abspath(tex_filename)
    latex_dir = os.path.dirname(abs_path)
    latex_filename = os.path.basename(abs_path)
    subprocess.call(["pdflatex", latex_filename], cwd=latex_dir)


def print_help():
    """
    Print the user help information to standard out.
    """
    print "USAGE: Call \"python peakCubeOverview.py HDF5File expIndex anaIndex outFile\" "
    print ""
    print "HDF5File : The OpenMSI HDF5 file to be used"
    print "expIndex : The index of the experiment within the file"
    print "anaIndex : The index of the analysis with the global peak finding results"
    print "outFile : The basename (and path) for the output files"
    print ""
    print "Example usage:"
    print "python test_peakcube.py /work2/bowen/2012_0403_KBL_platename_100212.h5 0 1 /work2/bowen/test_images/slice"
    print ""
    print "Output : "
    print "1) One image for each m/z slice of the peak_cube data"
    print "2) Latex file with a summary of all images. The latex file can be compiled to pdf using pdflatex."

if __name__ == "__main__":
    main()
