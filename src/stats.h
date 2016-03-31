#ifndef STATS_H
#define STATS_H

#include <vector>
#include <algorithm>

using namespace std;

double get_median(double * array, int length)
{
	vector<double> tmp ( array, array+length );
	sort( tmp.begin(), tmp.end() );
	return ( 0.5 * ( tmp[(length-1)/2] + tmp[length/2] ) );
}

double get_median_absolute_deviation(double * array, int length, double * median)
{
	*median = get_median( array, length );

	double * dev = new double [length];
	for (int i = 0; i < length; ++i)
		dev[i] = abs( array[i] - *median );

	double result = get_median( dev, length );

	delete [] dev;

	return ( result );
}

void find_outliers_Hampel(double * array, int length, bool * results, double * med)
{
	double mad = get_median_absolute_deviation(array, length, med);

	//use of this factor defines Hampel method (Hampel 1985)
	double factor = 5.2;

	//assumes results has already been allocated elsewhere
	for (int i = 0; i < length; ++i)
		results[i] = ( abs( array[i] - *med ) > factor * mad );	//is true if array[i] is an outlier (by this detection algorithm)

	return;
}

#endif

//End of file
