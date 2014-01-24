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


order = 1024
x,y = hilbert_curve(order=order)
val=np.arange(len(x))
img = np.zeros( shape=(order,order) )
img[x,y] = val
plt.imshow(img)
plt.show()
