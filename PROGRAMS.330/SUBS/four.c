void four(float data[], int nn, int isign, float *dt, float *df){
/*
	subroutine four(data,nn,isign,dt,df)
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df
c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(data(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  data is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     data are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array data, replacing the input data.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c-----
      real*4 data(*)
*/
	int n;
	int i, j, m, mmax, iiii, istep;
	float tempr, tempi;
	float wr, wi;
	float wstpr, wstpi;
	float sinth, theta;
	float dtt, dff;

	dtt = (*dt);
	dff = (*df);

	n = 2 * nn;
	if((dtt) == 0.0) dtt = 1./(nn*(dff)) ;
	if((dff) == 0.0) dff = 1./(nn*(dtt)) ;
	if((dtt) != (nn*(dff))) dff = 1./(nn*(dtt)) ;
	*dt = dtt;
	*df = dff;
	j = 1;
	for (i=1;i<= n ; i+=2){
		if(i < j){
			tempr = data[j-1];
			tempi = data[j  ];
			data[j-1] = data[i-1];
			data[j  ]=data[i  ];
			data[i-1] = tempr;
			data[i  ] = tempi;
		}
		m = n/2;
statement3:		if(j <= m) goto statement4 ;
			j = j-m;
			m = m/2;
			if(m >= 2)goto statement3 ;
statement4:		
		j=j+m;
    	}
	mmax = 2 ;
	while(mmax < n ){
		istep= 2 *mmax;
		theta = 6.283185307/(float)(isign*mmax);
		sinth=sin(theta/2.);
		wstpr=-2.*sinth*sinth;
		wstpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1; m <= mmax ; m +=2){
				for(i = m ; i <= n ; i+=istep){
					j=i+mmax;
					tempr=wr*data[j-1]-wi*data[j  ];
					tempi=wr*data[j  ]+wi*data[j-1];
					data[j-1]=data[i-1]-tempr;
					data[j  ]=data[i  ]-tempi;
					data[i-1]=data[i-1]+tempr;
					data[i  ]=data[i  ]+tempi;
				}
				tempr = wr;
				wr = wr*wstpr-wi*wstpi + wr;
				wi = wi*wstpr+tempr*wstpi + wi;
		}
		mmax = istep;
	}
 /*
  * 	get the correct dimensions for a Fourier Transform 
  * 	from the Discrete Fourier Transform
*/ 
	if(isign > 0){
		/*
		frequency to time domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*df);
		}
	} else {
		/*
		time to frequency domain
		*/
		for (iiii= 0 ; iiii < n ; iiii++){
			data[iiii] = data[iiii] * (*dt);
		}
 }
}
