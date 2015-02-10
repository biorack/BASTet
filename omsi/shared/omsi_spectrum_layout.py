"""This module provides capabilities for computing different layouts for spectra"""
import numpy as np
import math


def plot_2d_spectrum_as_image(hilbert_intensities, show_plot=False, show_axis=False):
    """
    Plot image with pixels colored according to hilbert_intensities.

    :param hilbert_intensities: 2D numpy array with the intensity values for the spectrum.
    :type hilbert_intensities: 2D numpy array.
    :param show_plot: Show the generated plot in a window.
    :type show_plot: Boolean
    :param show_axis: Show x,y axis for the plot. Default is False.
    :type show_axis: Boolean

    :returns: matplotlib image plot or None in case that the plotting failed.

    """
    try:
        if not show_plot:
            import matplotlib
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        if not show_plot:
            plt.ioff()
    except ImportError:
        return None
    # Plot the hilbert image
    plt.clf()
    if show_axis:
        imgplot = plt.imshow(np.log(hilbert_intensities + 1))
        imgplot.set_interpolation('nearest')
        if show_plot:
            plt.show()
        return imgplot, plt
    else:
        dpi = 50
        xsize = hilbert_intensities.shape[0] / float(dpi)
        ysize = hilbert_intensities.shape[1] / float(dpi)
        fig = plt.figure(frameon=False, figsize=(xsize, ysize), dpi=dpi)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(np.log(hilbert_intensities + 1),
                  aspect='normal', interpolation='nearest')
        if show_plot:
            plt.show()
        return ax, fig


def compute_hilbert_spectrum(original_coords, original_intensities, left=0, right=0):
    """
    Given a 1D spectrum, interpolate the spectrum onto the closest 2D hilbert curve.

    :param original_coords: The original coordinate values (m/z). Values must be increasing.
    :type original_coords: 1D numpy array in increasing order.
    :param original_intensities: The original intensity values. Same length as original_coords.
    :type original_intensities: 1D numpy array of same length as original_coords
    :param left: Optional. Value to be used for padding data at the lower bound during interpolation
    :type: float
    :param right: Optional. Value to be used for padding data at the upper bound during interpolation
    :type: float

    :returns: 2D numpy array with the coordinate (m/z) values for the hilbert spectrum and separate
             2D numpy array for the interpolated intensity values.

    :raises: ValueError If original_coords and original_intensities have different length.
    """
    # Check if the inputs are Ok
    if original_coords.shape[0] != original_intensities.shape[0]:
        raise ValueError("ERROR: omsi.shared.compute_hilbert_spectrum. The parameters original_coords and " +
                         "original_intensities must have same length.")

    # Compute the length of the two curves
    originallength = original_coords.shape[0]
    # hilbertlength =  math.pow(2, math.ceil(math.log(originallength)/math.log(2)))
    hilbertorder = math.sqrt(originallength)
    # find the closest power of 2
    hilbertorder = math.pow(2, math.ceil(math.log(hilbertorder) / math.log(2)))
    hilbertlength = hilbertorder * hilbertorder
    minmz = original_coords.min()
    maxmz = original_coords.max()
    rangemz = maxmz - minmz

    # Compute the hilbert curve
    hilbertmz = ((np.arange(hilbertlength, dtype=np.dtype('float'))
                  / float(hilbertlength - 1)) * rangemz) + minmz
    hilbertx, hilberty = hilbert_curve(order=hilbertorder)

    # Interpolate the intensity values
    hilbertintensities = reinterpolate_spectrum(
        coords=hilbertmz,
        original_coords=original_coords,
        original_intensitities=original_intensities,
        left=left,
        right=right)

    # Map the 1D hilbert mz and intensities values onto the correct 2D space
    outmz = np.zeros(shape=(hilbertorder, hilbertorder), dtype=hilbertmz.dtype)
    outmz[hilbertx, hilberty] = hilbertmz
    outintensities = np.zeros(
        shape=(hilbertorder, hilbertorder), dtype=hilbertmz.dtype)
    outintensities[hilbertx, hilberty] = hilbertintensities

    # Return the new hilbert spectrum
    return outmz, outintensities


def reinterpolate_spectrum(coords, original_coords, original_intensitities, left=0, right=0):
    """
    Given a 1D spectrum, interpolate the spectrum onto a new axis.

    :param coords: The coordinate values (m/z) for which intensities should be computed.
    :param original_coords: The original coordinate values (m/z). Values must be increasing.
    :param original_intensitities: The original intensity values. Same length as original_coords.
    :param left: Optional. Value to be used if coords < orignal_coords
    :param right: Optional. Value to be used if coords > orignal_coords

    :returns: y : {float, ndarray} The interpolated values, same shape as coords.

    :raises: ValueError If original_coords and original_intensities have different length.
    """
    return np.interp(x=coords,
                     xp=original_coords,
                     fp=original_intensitities,
                     left=left,
                     right=right)


def hilbert_curve(order=2):
    """
    Compute  a 2D hilbert curve.

    :param order: The order of the hilber curve. This is the length of the sides of the square, i.e.,
                 the number of points in x and y.
    :type order: Integer that defines a power of 2 (>=2)

    :returns: Returns two numpy arrays of integers x,y, indicating the locations of the
             vertices of the hilbert curve.

    """
    if order == 2:
        return np.array((0, 0, 1, 1)), np.array((0, 1, 1, 0))
    else:
        x, y = hilbert_curve(order / 2)
        xre = np.r_[y,
                    x,
                    order / 2 + x,
                    order - 1 - y]
        yre = np.r_[x,
                    order / 2 + y,
                    order / 2 + y,
                    order / 2 - 1 - x]
        return xre.astype('int'), yre.astype('int')
