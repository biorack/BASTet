import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array, exp, convolve
from collections import deque
import numpy as np

class findpeaks:
    Name = "findpeaks"
    def __init__(self, x, y, sizesmooth, slwindow, peakheight):
        self.x = x
        self.y = y
        self.smooth = sizesmooth
        self.slidewindow = slwindow
        self.peakheight = peakheight
        
    def smoothListGaussian(self):
        y=self.y
        size=self.smooth  

        x = np.array(range(-3*size,3*size))
        g = np.exp(-(x**2)/(2.0*float(size)**2))
        g = g / g.sum()
            
        smoothed = convolve(np.array(y),g, 'same')
        return smoothed
        
    def sliding_window_minimum(self):
        '''
        A iterator which takes the size of the window, `k`, and an iterable,
        `li`. Then returns an iterator such that the ith element yielded is equal
        to min(list(li)[max(i - k + 1, 0):i+1]).

        Each yield takes amortized O(1) time, and overall the generator takes O(k)
        space.
        '''
        __author__ = "Keegan Carruthers-Smith"
        __email__ = "keegan.csmith@gmail.com"
        __license__ = "MIT"
        li = self.y
        k = self.slidewindow
        window = deque()
        for i, x in enumerate(li):
            while window and window[-1][0] >= x:
                window.pop()
            window.append((x, i))
            while window[0][1] <= i - k:
                window.popleft()
            yield window[0][0]

    def peakdet(self):
        """
        Converted from MATLAB script at http://billauer.co.il/peakdet.html

        Currently returns two lists of tuples, but maybe arrays would be better

        function [maxtab, mintab]=peakdet(v, delta, x)
        %PEAKDET Detect peaks in a vector
        %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
        %        maxima and minima ("peaks") in the vector V.
        %        MAXTAB and MINTAB consists of two columns. Column 1
        %        contains indices in V, and column 2 the found values.
        %      
        %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
        %        in MAXTAB and MINTAB are replaced with the corresponding
        %        X-values.
        %
        %        A point is considered a maximum peak if it has the maximal
        %        value, and was preceded (to the left) by a value lower by
        %        DELTA.

        % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
        % This function is released to the public domain; Any use is allowed.

        """
        maxtab = []
        mintab = []
        v=self.y
        #x=self.x
        
        delta = self.peakheight
                
        #if x is None:
        x = arange(len(v))

        v = asarray(v)

        if len(v) != len(x):
            sys.exit('Input vectors v and x must have same length')

        if not isscalar(delta):
            sys.exit('Input argument delta must be a scalar')

        if delta <= 0:
            sys.exit('Input argument delta must be positive')

        mn, mx = Inf, -Inf
        mnpos, mxpos = NaN, NaN

        lookformax = True

        for i in arange(len(v)):
            this = v[i]
            if this > mx:
                mx = this
                mxpos = x[i]
            if this < mn:
                mn = this
                mnpos = x[i]

            if lookformax:
                if this < mx-delta:
                    maxtab.append((mxpos, mx))
                    mn = this
                    mnpos = x[i]
                    lookformax = False
            else:
                if this > mn+delta:
                    mintab.append((mnpos, mn))
                    mx = this
                    mxpos = x[i]
                    lookformax = True

        return maxtab, mintab

    def display(self):
        print self.Name
        print self.x
        print self.y
