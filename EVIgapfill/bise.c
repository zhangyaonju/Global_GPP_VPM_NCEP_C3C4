#include<stdio.h>
#include<stdlib.h>
//#include<R.h>
#include "bise.h"

//best index slope extraction
double *bise(int *rdays, double *ndvi, int *rslidperiod, double *corrected_ndvi)
{
	double slope_threshold = 0.2;
	double slope_thres_value = 0;
	double ndvi_chosen;
	double *temp_corrected_ndvi;
	double *temp_ndvi;
	int days,slidperiod;
	int i,j;
	int period, bypassed_elems;

	days = rdays[0];
	slidperiod = rslidperiod[0];

	temp_corrected_ndvi = (double *) calloc(3*days,sizeof(double));
	temp_ndvi = (double *) calloc(3*days,sizeof(double));
	
	//Triple data to avoid inconsistency of start and end of year
	for (i=0; i < days;i++ ){
		temp_ndvi[i] = ndvi[i];
		temp_ndvi[i+days] = ndvi[i];
		temp_ndvi[i+2*days] = ndvi[i];
	}

	//run algorithm
	temp_corrected_ndvi[0] = temp_ndvi[0];
	for ( i=1; i < (3*days); i++){
		if ( temp_ndvi[i] >= temp_ndvi[i-1] ){
			temp_corrected_ndvi[i] = temp_ndvi[i];
		} else {
			//Slope Extraction
			if ( (i+slidperiod) >= (3*days) )
				period = 3*days-i-1;
			else
				period = slidperiod;

			slope_thres_value = temp_ndvi[i] + slope_threshold * (temp_ndvi[i-1] - temp_ndvi[i]);

			bypassed_elems = 0;
			ndvi_chosen = 0;
			for ( j=i+1; j <= i+period-1; j++ ){
				//choose the next valid value within sliding period
				if (( temp_ndvi[j] > slope_thres_value ) && ( temp_ndvi[j] > ndvi_chosen )){
					ndvi_chosen = temp_ndvi[j];
					bypassed_elems = j-i; 
                		}
                		if (ndvi_chosen >= temp_ndvi[i-1])
                  			break;
			}
			if ( ndvi_chosen == 0 ) {
				temp_corrected_ndvi[i] = temp_ndvi[i];
			} else {
				//replace bypassed elements with -1
				for ( j=1; j <= bypassed_elems; j++)
					temp_corrected_ndvi[i-1+j] = -1;
				i = i + bypassed_elems;
				temp_corrected_ndvi[i] = ndvi_chosen;
			}
		}
	}
	
	//extract the middle of the 3 ndvi-copies
	for ( i=0; i<days; i++){
		corrected_ndvi[i] = temp_corrected_ndvi[i+days];
	}

	//free the allocated memory
	free(temp_corrected_ndvi);
	free(temp_ndvi);

	return corrected_ndvi;
}
