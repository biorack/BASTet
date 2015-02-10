import numpy as np
import matplotlib.pyplot as plt
 
def hilbert_curve(order=2):
    """Compute  a 2D hilbert curve.
 
       :param order: The order of the hilber curve. This is the length of the sides of the square, i.e., the number of points in x and y. 
       :type order: Integer that defines a power of 2 (>=2)
 
       :returns: Returns two arrays numpy arrays of integers x,y, indicating the locations of the vertices of the hilber curve. 
 
    """
    if order==2:
        return  np.array((0, 0, 1, 1)), np.array((0, 1, 1, 0))
    else:
        x, y = hilbert_curve(order/2)
        xre = np.r_[y, x,     order/2+x, order-1-y  ]
        yre = np.r_[x, order/2+y, order/2+y, order/2-1-x]
        return xre, yre

def plot_hilbert_image(x,y,val,order) : 
    """Plot image with pixels colored according to val and layed out 
       based on the hilbert coordinates given by x and y.
    """
    #Create hilbert image
    img = np.zeros( shape=(order,order) )
    img[x,y] = val

    #Plot the hilber image
    imgplot = plt.imshow(img)
    imgplot.set_interpolation('nearest')
    plt.show()

    return imgplot

def plot_hilbert_curve(x,y) :
    """Plot the hilbert curve specified by the coordinate arrays x,y."""
    plt.plot(x,y)
    plt.show()

def main(argv=None):
    """Then main function"""

    #Read input arguments
    import sys
    from sys import argv,exit
    if argv is None:
        argv = sys.argv
    order = 128
    if len(argv)>1 : 
        order = int(argv[1])

    #Compute the hilbert curve
    x,y = hilbert_curve(order=order)
    val=np.arange(len(x))
    
    plot_hilbert_image(x,y,val,order)
    plot_hilbert_curve(x,y)

if __name__ == "__main__":
    main()

