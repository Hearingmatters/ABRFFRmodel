#include <stdlib.h>
#include <string.h>
#include <math.h>
/*#include <minmax.h>*/

#include "mex.h"


/*  Generates event times from a Poisson process or a Poisson-equivalent integrate-and-fire process.
		
	Calling syntax:
		output = PoissonProc(rate, T)  for scalar 'rate' (homogeneous Poisson process)
		output = PoissonProc(rate, dt, method)  for vector 'rate' (inhomogeneous Poisson process)
	
	where		output		is a vector containing the times at which events occurred (sec).
				rate		is the rate vector or a scalar for a constant rate (events/sec).
				T			is the total simulation time (sec).
				dt			is the sample period of the rate function (sec).
				method		is a string specifying the simulation method [Default = 'approx']

		The possible method strings are:
			'approx'	 - Small time interval approximation method.  A point occurs at time t, where
							t is a multiple of dt, if U < rate(t) * dt, where U is a uniform random 
							variable and rate(t) is the value of the rate function at time t.
			'transform'	 - A unit rate homogeneous Poisson process is produced, and then time is
						    transformed via the rate function.
			'integ&fire' - Essentially the same as 'transform', which is equivalent to an integrate and
							fire process, except that negative "rate" values are set to zero in 'transform'.
							In the 'integ&fire' method, negative "rate" values are used instead of being
							set to zero.  Thus, if negative "rate" values are present, this method does
							not produce a Poisson process, although the resulting Poisson-like process
							is interesting and useful (see Jackson and Carney (2005) "The Spontaneous-Rate 
							Histogram of the Auditory Nerve Can Be Explained by Only Two or Three Spontaneous 
							Rates and Long-Range Dependence"  Journal of the Association for Research in
							Otolaryngology, 6(2), in press.
	
   Copyright � 2003-2005 by B. Scott Jackson
   Revision: 1.0    Date: March 14, 2005
   History:
       This file is a modification of a similar file 'poissproc.c' written by B. Scott Jackson.  The changes made
       in the creation of this file ('poissonproc.c') are:
       		- Additional documentation was added.
       		- The 'thinning' method was removed.
       		- The 'integ&fire' method was added and made equivalent to the previous 'transform' method, which
       			did not rectify (i.e. set negative values to zero) the input rate function.
       	    - The 'transform' method was modified to include rectification of the input rate function.
       	    - The 'approx' method was simplified so that simulation time step is always equal to the sampling period 
       	    	of the rate function.  Previously, a different simulation time step could be specified.

*/

enum method {None, IntervalRVs, DeltaTime, TimeTransform, IntegrateFire};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
mxArray		*randInputArray[1], *randOutputArray[1];
double		*randNums, *randDims;
long		randBufIndex;

double		*rate, T, dt;
double		meanRate, *output, time, *tempOutputPtr;
double		unitRateIntrvl, Xsum;
long		k, numRate, Nout, NoutMax, Nmore;
char		methodString[10];
enum method	methodType;

/* Check and get the input arguments. */
if ( (nrhs < 2) && (nrhs > 3) )  mexErrMsgTxt("Requires two or three input arguments.");

if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || ( (mxGetM(prhs[0]) != 1)
														&& (mxGetN(prhs[0]) != 1) ) )
	mexErrMsgTxt("The first argument must be a real, double-precision scalar or vector.");

rate = mxGetPr(prhs[0]);
numRate = mxGetNumberOfElements(prhs[0]);

if (numRate == 1)	/* If 'rate' is a scalar, then a homogeneous Poisson process will be simulated. */
	{
	if (nrhs > 2)
		mexErrMsgTxt("Exactly two input arguments required if 'rate' is a scalar.");
	
	if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetNumberOfElements(prhs[1]) != 1) )
		mexErrMsgTxt("The second argument must be a real, double-precision scalar.");
	
	methodType = IntervalRVs;
	T = mxGetScalar(prhs[1]);
	if (T <= 0)
		mexErrMsgTxt("The simulation time must be greater than zero.");
	
	dt = mxGetNaN(); /* 'dt' is not used for a homogeneous Poisson process. */
	}
else				/* If 'rate' is a vector, then an inhomogeneous Poisson process will be simulated. */
	{
	if ( (nrhs > 2) 
			&& ( !mxIsChar(prhs[2]) || mxIsEmpty(prhs[2]) || (mxGetNumberOfElements(prhs[2]) > 10) ) )
		mexErrMsgTxt("The third argument argument must be a valid method string.");

	if (nrhs < 3)
		methodType = DeltaTime;
	else
		{
		if ( mxGetString(prhs[2], methodString, 11) != 0 )
			mexErrMsgTxt("Could not read the method type string.");
		if (strcmp(methodString, "approx") == 0)
			methodType = DeltaTime;
		else if (strcmp(methodString, "transform") == 0)
			methodType = TimeTransform;
		else if (strcmp(methodString, "integ&fire") == 0)
			methodType = IntegrateFire;
		else
			mexErrMsgTxt("Invalid method type.");
		}
	
	if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetNumberOfElements(prhs[1]) != 1) )
		mexErrMsgTxt("The second argument must be a real, double-precision scalar.");
	
	dt = mxGetScalar(prhs[1]);
	if ( dt <= 0 )
		mexErrMsgTxt("The rate function sampling period must be greater than zero.");
	T = numRate * dt;
	}

/* Create the output buffer to store the spike times. */
meanRate = 0.0;
for (k=0; k<numRate; ++k)
	{
	meanRate += fmax(rate[k], 0.0);
	}
meanRate /= numRate;  /* Calculate the mean of the rate function. */

NoutMax = (long)fmax(ceil(meanRate * T), 100.0); /* Estimate the number of spikes that will occur. */
output = (double *)mxCalloc(NoutMax, sizeof(double));
Nout = 0;

/* Get a vector of pseudo-random numbers from MATLAB.  Calling MATLAB from a MEX-function is slow, 
	 so calling the MATLAB 'rand' function for each individual pseudo-random number needed greatly 
	 increases the computation time.																*/
randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
randDims = mxGetPr(randInputArray[0]);
randDims[0] = 1;
randDims[1] = (methodType == DeltaTime) ? numRate : NoutMax;
mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
randNums = mxGetPr(randOutputArray[0]);
randBufIndex = 0;

time = 0.0;

if (methodType == IntervalRVs)	/* Homogeneous Poisson process */
	{
	if (rate[0] > 0.0) /* If the "rate" is negative, an empty matrix will be returned. */
		{
		while (time < T)
			{
			if (Nout >= NoutMax)  /* If the output buffer for the spike times is too small . . . */
				{
				Nmore = (long)fmax(meanRate * (T-time), 100.0);
				NoutMax += Nmore;
				output = (double *)mxRealloc(output, NoutMax*sizeof(double));  /* . . . allocate additional memory . . . */
				if (output == NULL)  mexErrMsgTxt("Out of Memory");
				/* . . . and get more pseudo-random numbers since we will have used all of them up. */
				randDims[1] = Nmore;
				mxDestroyArray(randOutputArray[0]);
				mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
				randNums = mxGetPr(randOutputArray[0]);
				randBufIndex = 0;
				}
				
			time +=  (-log(randNums[randBufIndex++]))/rate[0];  /* Add the next interval length to the accumulated time . . . */
			output[Nout++] = time; 								/* . . . and output this time as the next spike time.         */
			}
		
		output[--Nout] = mxGetNaN(); /* Get rid of the last spike time, since it is after the end of the requested simulation time. */
		}
	}
else				/* Inhomogeneous Poisson process */
	{
	switch ( methodType )
		{
		case DeltaTime:
			/* Simulate the process. */
			for (k=0, time=0.0; k < numRate; ++k, time+=dt)
				{
				/* A spike occurs at this time with the probability: rate * dt. */
				if ( randNums[k] < rate[k] * dt )
					{
					output[Nout] = time;
					if (++Nout >= NoutMax) /* If the output buffer for the spike times is too small . . . */
						{
						/* . . . determine approximately how many more spikes we should expect, . . .*/
						NoutMax += (long)fmax(meanRate * (T-time), 100.0);
						output = (double *)mxRealloc(output, NoutMax*sizeof(double));  /* . . . and allocate additional memory. */
						if (output == NULL)  mexErrMsgTxt("Out of Memory");
						}
					}
				}
			break;
		
		case TimeTransform:
			/* Set negative "rate" values equal to zero. */
			for (k=0; k<numRate; ++k)
				if (rate[k] < 0)  rate[k] = 0;
		case IntegrateFire:
			/* Initialize process */
			unitRateIntrvl = -log(randNums[randBufIndex++])/dt;   /* An interval from a homogeneous Poisson process with a rate  */
																  /* of one divided by the sample period of the rate function.   */
																  /* This sample period is accounted for here, instead of in the */
																  /* integral (i.e. Xsum += rate[k] * dt), in order to reduce    */
																  /* the number of computations.                                 */
			Xsum = 0.0;
			
			/* Simulate the process.  This is mathematically equivalent to integrating the rate function (i.e. Xsum) until it is */
			/* larger than an interval from a homogeneous Poisson process (HPP) with a rate of one.  However, both sides of this */
			/* inequality are divided by the sampling period of the rate function, dt.                                           */
			for (k=0, time=dt; (k<numRate) && (time<T+dt); ++k, time+=dt)
				{
				Xsum += rate[k]; /* The integration process. */
				
				while (Xsum >= unitRateIntrvl) /* Once the integral is larger than the HPP interval, . . . */
					{
					Xsum -= unitRateIntrvl; /* ... reduce the integral by the length of the HPP interval and . . . */
					output[Nout] = time - Xsum*dt/rate[k];  /* . . . determine the time at which the process fires a spike, i.e. the */
															/*       exact time when the integral was equal to the HPP interval.     */
					if (++Nout >= NoutMax)  /* If the output buffer for the spike times is too small . . . */
						{
						/* . . . determine approximately how many more spikes we should expect, . . .*/
						Nmore = (long)fmax(meanRate * (T-output[Nout-1]), 100.0);
						NoutMax += Nmore;
						output = (double *)mxRealloc(output, NoutMax*sizeof(double));  /* . . . allocate additional memory, . . . */
						if (output == NULL)  mexErrMsgTxt("Out of Memory");
						/* . . . and get more pseudo-random numbers since we will have used all of them up. */
						randDims[1] = Nmore;
						mxDestroyArray(randOutputArray[0]);
						mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
						randNums = mxGetPr(randOutputArray[0]);
						randBufIndex = 0;
						}
					unitRateIntrvl = -log(randNums[randBufIndex++])/dt;  /* Get another interval from an HPP and divide it by the */
																	     /* sample period of the rate function (see above).       */
					}
				}
			
			for (; (Nout>0)&&(output[Nout-1]>T); --Nout)  output[Nout-1] = mxGetNaN(); /* Remove spike times later than the requested */
																					   /* simulation end time.                        */ 
			break;
		}
	}


mxDestroyArray(randInputArray[0]);
mxDestroyArray(randOutputArray[0]);

/* Finish up */
if (Nout > 0)
	{
	/* Deallocate the unused portion of the output array. */
	tempOutputPtr = (double *)mxRealloc(output, Nout*sizeof(double));
	if (tempOutputPtr == NULL)
		{
		mexPrintf("Nout = %d; NoutMax = %d\n", Nout, NoutMax);
		mexWarnMsgTxt("Unable to deallocate the unused portion of memory in the output array.");
		for (k=Nout; k<NoutMax; ++k)
			output[k] = mxGetNaN();
		Nout = NoutMax;
		}
	else
		output = tempOutputPtr;
	
	/* Create the output and set its pointer to the output array. */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxFree(mxGetPr(plhs[0]));
	mxSetN(plhs[0], Nout);
	mxSetPr(plhs[0], output);
	}
else  /* Return an empty matrix. */
	{
	mxFree(output);
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
	}
}
