
// ----

int my2Dindex(int xy, int z, int sizeZ) {
 	// return z + (sizeZ * y) + (sizeZ * sizeY * x);
	return z + (sizeZ * xy);
}

__kernel void clSmoothed( __global const clfloat *input,
						__global const clfloat *g,
						int smoothsize,
						int xysize, int zsize,
                  		__global clfloat *output
						)
{
	int gid = get_global_id(0);
	// int gid2 = get_global_id(1);
    
	int myindex;
	int i;
	
	int a_range = 3*smoothsize;
	int a_size = a_range*2;
	int c_size = a_size + zsize - 1;
	int kmin, kmax, k;
	
	clfloat tempvar;
	int olimiter = 0;
	
	// convolve
	for (i = 0; i < c_size; i++)
	{
		tempvar = 0;
		if(i >= (zsize - 1)) {
			kmin = i - (zsize - 1);
		} else {
			kmin = 0;
		}
		
		if( i < (a_size - 1)) {
			kmax = i;
		} else {
			kmax = a_size - 1;
		}
		
		for (k = kmin; k <= kmax; k++)
		{
			myindex = my2Dindex(gid,i-k,zsize);
			tempvar = tempvar + g[k] * input[myindex];
		}

		if(i >= (a_range-1) && olimiter < zsize) {
			myindex = my2Dindex(gid,i - (a_range-1),zsize);
			output[myindex] = tempvar;
			olimiter++;
		}
	}
}

__kernel void clSliding_window_minimum(	__global const clfloat *li,
									__global clfloat *windowValue, __global int *windowIndex,
									int wsize,
									int k,
									int xysize, int zsize,
									__global clfloat *slmin
									)
{
    int gid = get_global_id(0);

	int i = 0;
	int myindex;
		
	int empty = 0;
	int firstindex = 0;
	clfloat firstvalue = 0;
	clfloat lastvalue = 0;
	clfloat livalue = 0;
	int count = 0;

	int front=0;
	int rear=-1;
	
	for(i=0; i < zsize; i++) {
		// empty = isempty();
		myindex = my2Dindex(gid,i,zsize);
		livalue = li[myindex];
		
		if( count == 0 ) {
			empty = 1;
		} else {
			empty = 0;
			myindex = my2Dindex(gid,rear,wsize);
			lastvalue = windowValue[myindex];		// lastvalue = lastelementValue();
			// lastvalue = window[rear].Value;		// lastvalue = lastelementValue();
			
		}
				
		while(empty == 0 && lastvalue >= livalue) {
			// deleterear();
			rear = rear - 1;
			count--;
			
			// empty = isempty();
			if( count == 0 ) {
				empty = 1;
			} else {
				empty = 0;
				myindex = my2Dindex(gid,rear,wsize);
				lastvalue = windowValue[myindex];			// lastvalue = lastelementValue();
				// lastvalue = window[rear].Value;		// lastvalue = lastelementValue();
			}
		}
		
		// insertrear(li[i],i);
		rear = rear + 1;
		
		myindex = my2Dindex(gid,rear,wsize);
		windowValue[myindex] = livalue;
		windowIndex[myindex] = i;
		// window[rear].Value = livalue;
		// window[rear].Index = i;
		count++;
	
		// firstindex = firstelementIndex();
		myindex = my2Dindex(gid,front,wsize);
		firstindex = windowIndex[myindex];
		// firstindex = window[front].Index;
		
		while( firstindex <= (i - k) ) {
			// deletefront();
			front = front + 1;
			count--;
			myindex = my2Dindex(gid,front,wsize);
			firstindex = windowIndex[myindex];
			// firstindex = window[front].Index;
		}
		
		// firstvalue = firstelementValue();
		myindex = my2Dindex(gid,front,wsize);
		firstvalue = windowValue[myindex];
		// firstvalue = window[front].Value;
		myindex = my2Dindex(gid,i,zsize);
		slmin[myindex] = firstvalue;
	}
}


__kernel void clPeakdet( __global const clfloat *v,
	__global clfloat *my_maxtabVal,	__global int *my_maxtabPos, __global int *my_maxtabTot,
	__global clfloat *my_mintabVal, __global int *my_mintabPos, __global int *my_mintabTot,
	int delta,
	int xysize, int zsize,
	int my_maxt_size, int my_mint_size
	)
{

	int gid = get_global_id(0);

	int i, myindex;

	// initial values (first z value at an x/y pos)
	myindex = my2Dindex(gid,0,zsize);
	
	clfloat mn = v[myindex];
	clfloat mx = v[myindex];
	clfloat my_v = 0;	
	int mnpos = 0;
	int mxpos = 0;
	int mincounter = 0;
	int maxcounter = 0;
	int lookformax = 1;			// True (1), False(0) : Flag	
	
	for(i = 0; i < zsize; i++) {
		myindex = my2Dindex(gid,i,zsize);
		my_v = v[myindex];
		
		if(my_v > mx) {
			mx = my_v;
			mxpos = i;
		}
		
		if(my_v < mn) {
			mn = my_v;
			mnpos = i;
		}
		
		if(lookformax == 1) {
			if(my_v < (mx - delta)) {
				myindex = my2Dindex(gid,maxcounter,my_maxt_size);
				my_maxtabVal[myindex] = mx;
				my_maxtabPos[myindex] = mxpos;
				maxcounter++;
				mn = my_v;
				mnpos = i;
				lookformax = 0;
			}
		} else {
			if (my_v > (mn + delta)) {
				myindex = my2Dindex(gid,mincounter,my_mint_size);
				my_mintabVal[myindex] = mn;
				my_mintabPos[myindex] = mnpos;
				mincounter++;
				mx = my_v;
				mxpos = i;
				lookformax = 1;
			}
		}
		
	}
	
	myindex = my2Dindex(gid,0,1);
	my_maxtabTot[myindex] = maxcounter;
	my_mintabTot[myindex] = mincounter;
}

