"""This module provides functionality for reading img mass spectrometry image files"""
import os
import numpy as np
import warnings
from omsi.dataformat.file_reader_base import file_reader_base


class img_file(file_reader_base):

    """Interface for reading a single 2D img file

       The img format consists of three different files:
       i) hdr header file, ii) t2m which contains the m/z data,
       iii) img data file.
    """

    def __init__(self, hdr_filename=None, t2m_filename=None, img_filename=None, basename=None, readdata=True):
        """Open an img file for data reading.

            :param hdr_filename: The name of the hdr header file
            :type hdr_filename: string

            :param t2m_filename: The name of the t2m_filename
            :type t2m_filename: string

            :param img_filename: The name of the img data file
            :type img_filename: string

            :param basename: Instead of img_filename, t2m_filename, and hdr_filename one may also supply just
                             a single basename. The basename is completed with the .img, .t2m, .hdr extension
                             to load the data.
            :type basename: string

            :param readdata: Unused here. Included as this is a required parameter for the base API.
            :type readdata: Boolean

            :raises ValueError: In case that basename and hdr_filename, t2m_filename, and img_filename are specified.
        """
        super(img_file, self).__init__(basename, readdata)
        self.data_type = 'uint16'
        self.shape = [0, 0, 0]  # Number of pixels in x,y, and z. NOTE: Type changed to tuple later on.
        self.mz = 0  # A numpy vector with the m/z values of the instrument

        if basename and hdr_filename and t2m_filename and img_filename:
            raise ValueError(
                "Conflicting input. Provide either basename or the " +
                "hdr_filename,t2m_filename,img_filename parameters but not both.")
        if basename:
            basefile = basename
            if os.path.isdir(basename):
                filelist = self.get_files_from_dir(basename)
                print filelist
                if len(filelist) > 0:
                    basefile = filelist[0]
                else:
                    raise ValueError("No valid img file found in the given directory.")
            print basefile
            if os.path.exists(basefile + ".hdr") and \
                    os.path.exists(basefile + ".t2m") and \
                    os.path.exists(basefile + ".img"):
                hdr_filename = basefile + ".hdr"
                t2m_filename = basefile + ".t2m"
                img_filename = basefile + ".img"
            else:
                raise ValueError("No valid img file found for the given basename.")
        elif hdr_filename and t2m_filename and img_filename:
            pass  # Nothing to be done
        else:
            raise ValueError("Missing input parameter. Either provide: " +
                             " i) basename or ii) hdr_filename, t2m_filename, img_filename")
        
        # Initialize the x and y length
        hdr = open(hdr_filename, 'rb')
        hdrdata = np.fromfile(file=hdr_filename, dtype='int16', count=-1)
        self.shape[0] = hdrdata[23]
        self.shape[1] = hdrdata[22]
        hdr.close()

        # Initialize the z length
        t2m = open(t2m_filename, 'rb')
        self.mz = np.fromfile(file=t2m, dtype='float32', count=-1)
        self.shape[2] = self.mz.shape[0]
        t2m.close()

        #Convert the shape variable to the expected tuple
        self.shape = tuple(self.shape)

        # Open the img file with the spectrum data
        self.img_filename = img_filename
        self.file_opened = False
        try:
            self.m_img_file = np.memmap(filename=self.img_filename,
                                        dtype=self.data_type,
                                        shape=self.shape,
                                        mode='r',
                                        order='C')
            self.file_opened = True
        except ValueError:
            # Check if the size of the file matches what we expect
            imgsize = os.stat(self.img_filename).st_size
            itemsize = np.dtype(self.data_type).itemsize
            expectednumvalues = int(self.shape[0]) * int(self.shape[1]) * int(self.shape[2]) 
            expectedsize = expectednumvalues * int(itemsize)
            sizedifference = expectedsize - imgsize
            print "IMG size: " + str(imgsize) + " Expected size: " + \
                  str(expectedsize) + "  (difference="+str(sizedifference) + ")"
            if imgsize < expectedsize:
                # Check whether the missing data aligns with images or spectra
                slicesize = int(self.shape[0]) * int(self.shape[1]) * itemsize
                spectrumsize = int(self.shape[2]) * itemsize
                percentmissing = float(sizedifference)/float(expectedsize)
                valuesmissing = float(sizedifference) / itemsize
                warnings.warn("WARNING: Missing "+str(sizedifference) +
                              " bytes in img file (missing " + str(valuesmissing) +
                              " intensity values; "+str(percentmissing)+"%)." + 
                              " Expected shape: "+str(self.shape)) 
                # Define how we should deal with the error
                expandslice = (sizedifference % slicesize) == 0
                expandspectra = (sizedifference % spectrumsize) == 0
                if not expandslice:
                    expandspectra = True
                # Complete missing spectra
                if expandspectra:
                    warnings.warn("Dealing with missing data in img file by completing last spectra with 0's.")
                    # TODO np.require create an in-memory copy of the full data. Provide option to use memmap'ed tempfile.
                    tempmap = np.require(np.memmap(filename=self.img_filename,
                                                   dtype=self.data_type,
                                                   mode='r',
                                                   order='C'),
                                         requirements=['O', 'C'])
                    #Extend the memmap to the expected size
                    tempmap.resize((expectednumvalues, ))
                    # Reshape the memmap to the expected shape
                    self.m_img_file = tempmap.reshape(self.shape, order='C')
                    self.file_opened = True
                # Complete missing slices
                elif expandslice:
                    slicesmissing = sizedifference / slicesize
                    self.mz = self.mz[:(-slicesmissing)]
                    warnings.warn("Dealing with missing data in img file by updating he m/z axis.." +
                                  " It looks like the m/z axis data may be inconsistent" +
                                  " with the binary data. Removing "+str(slicesmissing) +
                                  " bins from the m/z axis.")
                    self.shape = list(self.shape)
                    self.shape[2] = self.mz.shape[0]
                    self.shape = tuple(self.shape)
                    self.m_img_file = np.memmap(filename=self.img_filename,
                                                dtype=self.data_type,
                                                shape=self.shape,
                                                mode='r',
                                                order='C')
                    self.file_opened = True
                else:
                    raise
            else:
                raise
        except:
            print "Error while opening the img file: " + img_filename
            raise

    def __getitem__(self, key):
        """Enable slicing of img files"""
        if self.m_img_file is None:
            self.m_img_file = np.memmap(
                filename=self.img_filename, dtype=self.data_type, shape=self.shape, mode='r', order='C')
        return self.m_img_file[key]

    def close_file(self):
        """Close the img file"""
        if self.file_opened:
            del self.m_img_file
            self.m_img_file = None
            self.file_opened = False

    @classmethod
    def is_valid_dataset(cls, name):
        """Check whether the given file or directory points to a img file.

           :param name: Name of the dir or file.
           :type name: unicode

           :returns: Boolean indicating whether the given file or folder is a valid img file.
        """
        # Check if this points to a basname for the img file
        if os.path.exists(name + ".hdr") and os.path.exists(name + ".t2m") and os.path.exists(name + ".img"):
            return True
        # If we point to a director, check if the dir contains an img file
        elif os.path.isdir(name):
            filelist = cls.get_files_from_dir(name)
            if len(filelist) > 0:
                return True

        return False

    @classmethod
    def size(cls, name):
        """
        Classmethod used to check the estimated size for the given file/folder.

        :param name: Name of the dir or file.
        :type name: unicode

        :returns: Integer indicating the size in byte or None if unknown.
        """
        basename = None
        if os.path.exists(name + ".hdr") and os.path.exists(name + ".t2m") and os.path.exists(name + ".img"):
            basename = name
        # If we point to a director, check if the dir contains an img file
        elif os.path.isdir(name):
            filelist = cls.get_files_from_dir(name)
            if len(filelist) > 0:
                basename = filelist[0]
        if basename is not None:
            return os.stat(name + ".img").st_size
        else:
            return None

    @classmethod
    def get_files_from_dir(cls, dirname):
        """
        Get a list of all basenames of all img files in a given directory.
        Note: The basenames include the dirname.
        """
        filelist = []
        for l in os.listdir(dirname):
            currname = os.path.join(dirname, l)
            if os.path.isfile(currname) and \
                    currname.endswith(".img"):
                basename = currname.rstrip(".img")
                if os.path.exists(basename + ".hdr") and \
                        os.path.exists(basename + ".t2m") and \
                        os.path.exists(basename + ".img"):
                    filelist.append(basename)
        return filelist

    def __del__(self):
        """Close the file before garbage collection"""
        self.close_file()
