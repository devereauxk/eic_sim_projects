#include <iostream>

void logbins(double xmin, double xmax, int nbins, double xbins[])
{
	double xmin_log = log10(xmin);
	double xmax_log = log10(xmax);
	double xbinw_log = (xmax_log-xmin_log)/nbins; // bin width in log scale
	xbins[0] = xmin;
	for (int ibin = 1; ibin < nbins; ++ibin)
	{
		xbins[ibin] = pow(10, xmin_log+ibin*xbinw_log);
	}
	xbins[nbins] = xmax;

	cout << "Log bins succesfully initialized: {";
	for (int ibin = 0; ibin < nbins; ++ibin)
	{
		cout << xbins[ibin] << ", ";
	}
	cout << xbins[nbins] << "};" << endl;
}

void linbins(double xmin, double xmax, int nbins, double xbins[])
{
	double xbinw = (xmax-xmin)/nbins; // bin width in linear scale
	xbins[0] = xmin;
	for (int ibin = 1; ibin < nbins; ++ibin)
	{
		xbins[ibin] = xmin+ibin*xbinw;
	}
	xbins[nbins] = xmax;

	cout << "Linear bins succesfully initialized: {";
	for (int ibin = 0; ibin < nbins; ++ibin)
	{
		cout << xbins[ibin] << ", ";
	}
	cout << xbins[nbins] << "};" << endl;
}