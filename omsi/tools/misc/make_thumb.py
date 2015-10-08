"""Simple script to generate thumbnail images"""

from omsi.dataformat.omsi_file.main_file import omsi_file

try:
    from PIL import Image
except ImportError:
    import Image
import numpy as np
import sys


def main(argv=None):
    """Then main function"""
    if argv is None:
        argv = sys.argv

    # Check for correct usage
    if len(argv) < 5:
        print "Usage: python makeThumb.py <HDF5-File> <NMF-Index1> <NMF-Index2> <NMF-Index3>"
        sys.exit(0)

    # Settings
    nmf_analysis_index = 1
    experiment_index = 0
    input_filename = argv[1]
    image_index1 = int(argv[2])
    image_index2 = int(argv[3])
    image_index3 = int(argv[4])
    apply_log_scale = True
    thumbnail_filename = input_filename + "_0.png"

    # Open the file and required analysis dataset
    omsi_input_file = omsi_file(input_filename, 'r')
    input_experiment = omsi_input_file.get_experiment(experiment_index)
    input_analysis = input_experiment.get_analysis(nmf_analysis_index)
    ho_data = input_analysis['ho']

    # Load the required data
    num_pixel_x = ho_data.shape[0]
    num_pixel_y = ho_data.shape[1]
    image_data1 = ho_data[:, :, image_index1].reshape((num_pixel_x, num_pixel_y))
    image_data2 = ho_data[:, :, image_index2].reshape((num_pixel_x, num_pixel_y))
    image_data3 = ho_data[:, :, image_index3].reshape((num_pixel_x, num_pixel_y))

    # Scale the data
    if apply_log_scale:
        image_data1 = np.log(image_data1 + 1)
        image_data2 = np.log(image_data2 + 1)
        image_data3 = np.log(image_data3 + 1)

    # Normalize the data values
    image_data1 /= float(np.max(image_data1))
    image_data2 /= float(np.max(image_data2))
    image_data3 /= float(np.max(image_data3))

    # Generate the grayscale images
    image_channel_red = Image.fromarray(image_data1.astype('float') * 255).convert('L')
    image_channel_green = Image.fromarray(image_data2.astype('float') * 255).convert('L')
    image_channel_blue = Image.fromarray(image_data3.astype('float') * 255).convert('L')

    # Generate the RGB image and save the file
    thumbnail = Image.merge('RGB', (image_channel_red, image_channel_green, image_channel_blue))
    print "Save image"
    thumbnail.save(thumbnail_filename, 'PNG')
    print thumbnail_filename


if __name__ == "__main__":
    main()

