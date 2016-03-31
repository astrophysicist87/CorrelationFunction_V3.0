#ifndef STATS_H
#define STATS_H

#include <vector>
#include <algorithm>

using namespace std;

double get_median(double * array, int length);
double get_median_absolute_deviation(double * array, int length, double * median);
void find_outliers_Hampel(double * array, int length, bool * results, double * med);

#endif

//End of file
