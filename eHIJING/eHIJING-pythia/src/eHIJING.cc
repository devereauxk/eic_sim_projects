#include "eHIJING/eHIJING.h"
#include <iostream>
#include <fstream>
#include "eHIJING/integrator.h"
#include <cmath>
#include <thread>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#include <string>
#include <sstream>
#include <filesystem>
#include <atomic>

namespace EHIJING {

// color algebra constnats
const double CA = 3;
const double dA = 8;
const double CF = 4./3.;
// Lambda QCD and Lambda^2
const double mu = 0.25; // [GeV]
const double mu2 = std::pow(mu, 2); // GeV^2
// minimum and maximum nuclear thickness function
// for tabulation
const double TAmin = 0.1/5.068/5.068;
const double TAmax = 2.8/5.068/5.068;
// The b parameter for 3-flavor QCD
// b0 = 11 - 2/3*Nf
const double b0 = 9./2.;
const double twoPioverb0 = 2.*M_PI/b0;
const double piCAoverdA = M_PI * CA / dA;
const double Mproton = 0.938;
const double rho0 = 0.17/std::pow(5.068,3);
const double r0 = 1.12*5.068;
// coupling in eHIJING
double alphas(double Q2){
    double Q2overmu2 = Q2/mu2;
    if (Q2overmu2<2.71828) return twoPioverb0;
    else return twoPioverb0/std::log(Q2/mu2);
}

double alphas_PhiG(double x, double q2, double Qs2,
		   double powerG, double lambdaG, double K){
    return K*std::pow(1.-x, powerG)*std::pow(x, lambdaG) / (Qs2+q2);
}

// integrate du (1-cos(1/u)) / u
double inte_C(double u){
    return gsl_sf_Ci(1./u) + std::log(u);
}

double CHT_F1(double x){
    return 1.0 - sin(x)/x;
}

///////// Class: Multiple Collision /////////////////////
// initializer
MultipleCollision::MultipleCollision(double Kfactor, double powerG, double lambdaG):
Kfactor_(Kfactor),
powerG_(powerG),
lambdaG_(lambdaG),
rd(),
gen(rd()),
flat_gen(0.,1.),
Qs2Table(3, {31,31,31}, // TA, ln(x), ln(Q2) --> size and grid for Qs table
           {TAmin, std::log(1e-6), std::log(1.0)},
           {TAmax, std::log(.99), std::log(1e5)}
       )
{
}
// Tabulate Qs
void MultipleCollision::Tabulate(std::filesystem::path table_path){
    std::filesystem::path fname = table_path/std::filesystem::path("Qs.dat");
    if (std::filesystem::exists(fname)) {
        std::cout << "Loading Qs Table" << std::endl;
        std::ifstream f(fname.c_str());
        int count = 0;
        double entry;
        std::string line;
        while(getline(f, line)){
            std::istringstream in(line);
            in >> entry;
            if (count>=Qs2Table.size()){
                std::cerr << "Loading table Qs: mismatched size - 1" << std::endl;
                exit(-1);
            }
            Qs2Table.set_with_linear_index(count, entry);
	    count ++;
        }
        if (count<Qs2Table.size()){
            std::cerr << "Loading table Qs: mismatched size - 2" << std::endl;
            exit(-1);
        }
    }
    else {
        if (std::filesystem::create_directory(table_path) ) std::cout << "Generating Qs^2 table" << std::endl;
        else std::cout << "Create dir failed" << std::endl;
        std::ofstream f(fname.c_str());

        // Table Qs as a function of lnx, lnQ2, TA
        for (int c=0; c<Qs2Table.size(); c++) {
            auto index = Qs2Table.LinearIndex2ArrayIndex(c);
            auto xvals = Qs2Table.ArrayIndex2Xvalues(index);
            double TA = xvals[0];
            double xB = std::exp(xvals[1]);
	    double Q2 = std::exp(xvals[2]);
            double aQs2 = compute_Qs2(TA, xB, Q2);
            Qs2Table.set_with_linear_index(c, aQs2);
            f << aQs2 << std::endl;
        }
    }
}
// self-consisten euqation for Qs2
double MultipleCollision::Qs2_self_consistent_eq(double Qs2, double TA, double xB, double Q2){
    double LHS = 0.;
    double scaledTA = piCAoverdA * TA;
    auto dfdq2 = [this, Qs2, xB, Q2](double ln1_q2oQs2) {
        double q2 = Qs2*(std::exp(ln1_q2oQs2)-1);
        double Jacobian = Qs2+q2;
        double xg = q2/Q2*xB;
        return alphas_PhiG(xg, q2, Qs2, this->powerG_, this->lambdaG_, this->Kfactor_) * Jacobian;
    };
    double error;
    double res =  scaledTA * quad_1d(dfdq2, {std::log(1.+.01*mu2/Qs2), std::log(1.+Q2/xB/Qs2)}, error);
    return res - Qs2;
}
// Solver of the Qs2 self-consistent equation, using a simple bisection
double MultipleCollision::compute_Qs2(double TA, double xB, double Q2){
    // a naive bisection
    double xleft = mu2, xright = Q2*100;
    const double EPS = 1e-4;
    double yleft = Qs2_self_consistent_eq(xleft, TA, xB, Q2),
           yright = Qs2_self_consistent_eq(xright, TA, xB, Q2);
    if (yleft*yright>0) {
        std::cout << "eHIJING warning: setting Qs2 = mu2" << std::endl;
        return xleft;
    } else {
        do {
            double xmid = (xright+xleft)/2.;
            double ymid = Qs2_self_consistent_eq(xmid, TA, xB, Q2);
            if (yleft*ymid<0) {
                yright = ymid;
                xright = xmid;
            } else{
                yleft = ymid;
                xleft = xmid;
            }
        } while (xright-xleft > EPS);
        return (xright+xleft)/2.;
    }
}

// Sample all elastic collisions, without radiation, ordered from high scale to low
double RealisticDensityProfile(double r, int A){
    double RA = EHIJING::r0*std::pow(A, 1./3.);
    double a0 = 0.5*5.068;
    // Woods Saxon profile
    double F = 1./(1.+std::exp((r-RA)/a0));
    double F0 = 1./(1.+std::exp(-RA/a0));
    // Normalize the profile to unity at r=0 (tiny effect when RA>>a0).
    return F/F0;
}

int MultipleCollision::sample_all_qt2(int pid, double E, double L, double Thickness, double xB, double Q2,
                             std::vector<double> GeoInfo,
		                     std::vector<double> & q2_list, std::vector<double> & t_list,
		                     std::vector<double> & phi_list) {
    q2_list.clear();
    t_list.clear();
    phi_list.clear();
    double TA = rho0*L;
    double CR = (pid==21)? CA : CF;
    double qs2 = TA / Thickness * Qs2(Thickness, xB, Q2);
    double tildeTA = Kfactor_*M_PI*CR*TA/dA;
    if (tildeTA<1e-6) return q2_list.size();
    double qs2overCTA = qs2/tildeTA;
    double q2max = Q2/xB;
    double q2min = .01*EHIJING::mu2;
    double q2 = q2max;
    double t0 = EHIJING::r0; // exlucde double scattering on the same nucleon
    double R0sq = GeoInfo[0]; 
    double V2 = GeoInfo[1];
    double TwoR0dotV = GeoInfo[2];
    double A = GeoInfo[3];
    if (L<t0) return q2_list.size();
    while (q2 > q2min) {
        // sample the next hard multiple collision
        double lnr = std::log(flat_gen(gen));
        double xg = q2/q2max;
	q2 = q2 * std::pow(1.0 + (lambdaG_-1.)*lnr*q2/tildeTA*std::pow(xg, -lambdaG_),
                                   1./(lambdaG_-1.));
        double t = t0 + flat_gen(gen)*(L-t0);
        double Rdistance = std::sqrt(R0sq + V2*t*t + TwoR0dotV*t);
	double DensityProfielRejection = RealisticDensityProfile(Rdistance, A);
        if ( flat_gen(gen) < std::pow(1.-xg, powerG_) * (q2/(q2+qs2)) * DensityProfielRejection ) {
	    // correct for the gluon distribtuion at large xg,
	    // and the screening effect
            // and any deviation from a constant nucleon density
            q2_list.push_back(q2);
            t_list.push_back(t);
            phi_list.push_back(2.*M_PI*flat_gen(gen));
	}
    }
    return q2_list.size();
}

/////////// Class: eHIJING
// initializer
eHIJING::eHIJING(int mode, double Kfactor,
                 double powerG, double lambdaG):
MultipleCollision(Kfactor, powerG, lambdaG),
mode_(mode),
Kfactor_(Kfactor),
powerG_(powerG),
lambdaG_(lambdaG),
rd(),
gen(rd()),
flat_gen(0.,1.),
GHT_Angular_Table(2, {51, 51}, // X = delta = 2kq/(k^2+q^22),
                               // ln(1+Y) = log(1+ (|k|-|q|)^2*t/(2*z*(1-z)*E) )
                               // 0.5<delta<1, ln(1)<ln(1+Y)<ln(11)
           {0.5, 1e-3},
           {0.99, std::log(11.0)}
       )
{
}
// Tabulate the Qs and (if necessary) the generlized H-T / GLV table (collinear H-T do not need separate table)
void eHIJING::Tabulate(std::filesystem::path table_path){
    MultipleCollision::Tabulate(table_path);
    if (mode_==1) {
        // Generlized higher-twist able is at least 4-dimensional with each entry a 2D integral
        // the following routine will use the max number of hard ware concurrency of your computer
        // to parallel the computation of the table
        std::filesystem::path fname = table_path/std::filesystem::path("GHT.dat");
        if (std::filesystem::exists(fname)) {
            std::cout << "Loading GHT Table" << std::endl;
            std::ifstream f(fname.c_str());
            int count = 0;
            double entry1, entry2;
            std::string line;
            while(getline(f, line)){
                std::istringstream in(line);
                in >> entry1;
                if (count>=GHT_Angular_Table.size()){
                    std::cerr << "Loading table GHT: mismatched size - 1" << std::endl;
                    exit(-1);
                }
                GHT_Angular_Table.set_with_linear_index(count, entry1);
                count ++;
            }
            if (count<GHT_Angular_Table.size()){
                std::cerr << "Loading table GHT: mismatched size - 2" << std::endl;
                exit(-1);
            }
        } else {
            std::filesystem::create_directory(table_path);
            std::ofstream f(fname.c_str());
            // Table Qs as a function of lnx, lnQ2, TA
            std::atomic_int counter =  0;
            int percentbatch = int(GHT_Angular_Table.size()/100.);
            auto code = [this, percentbatch](int start, int end) {
                static std::atomic_int counter;
                for (int c=start; c<end; c++) {
                    counter ++;

                    if (counter%percentbatch==0) {
                      std::cout <<std::flush << "\r" << counter/percentbatch << "% done";
                    }
                    auto index = GHT_Angular_Table.LinearIndex2ArrayIndex(c);
                    auto xvals = GHT_Angular_Table.ArrayIndex2Xvalues(index);
                    double X = xvals[0];
                    double ln1Y = xvals[1];
                    double Y = std::exp(ln1Y)-1.;
                    double entry1 = 0.;
                    entry1 = compute_GHT_Angular_Table(X, Y);
                    GHT_Angular_Table.set_with_linear_index(c, entry1);
                }
            };
            std::vector<std::thread> threads;
            int nthreads = std::thread::hardware_concurrency();
            int padding = int(std::ceil(GHT_Angular_Table.size()*1./nthreads));
            std::cout << "Generating GHT angular tables with " << nthreads << " thread" << std::endl;
            for(auto i=0; i<nthreads; ++i) {
                int start = i*padding;
                int end = std::min(padding*(i+1), GHT_Angular_Table.size());
                threads.push_back( std::thread(code, start, end) );
            }
            for(auto& t : threads) t.join();
            for (int c=0; c<GHT_Angular_Table.size(); c++) {
                f << GHT_Angular_Table.get_with_linear_index(c) << std::endl;
            }
        }
        std::cout << "... done" << std::endl;
    }
}

// computation of generalized HT table
double eHIJING::compute_GHT_Angular_Table(double X, double Y) {
    double A = Y/(1.-X);
    double B = Y*X/(1.-X);
    auto dfdphi = [X, Y, A, B](double phi) {
        double cosphi = std::cos(phi);
        double xcphi = X*cosphi;
        return xcphi/(1.-xcphi) * (1.-std::cos(A-B*cosphi));
    };
    double error;
    double result = quad_1d(dfdphi, {0., M_PI}, error);
    result /= M_PI;
    return result;
}

bool eHIJING::next_kt2_stochastic(double & kt2, 
	                int pid,
                        double E,
                        double kt2min,
                        std::vector<double> qt2s,
                        std::vector<double> ts) {
    double CR = (pid==21)? CA : CF;
    double CAoverCR = CA/CR;
    double CR_2overb0 = CR*2.0/b0;
    double zmin = std::min(.4/E, .4);
    double zmax = 1. - zmin;
    double logvac = std::log(zmax/zmin);
    int Ncolls = ts.size();
    double acceptance = 0.;
    if (mode_ == 0){
        while (acceptance<flat_gen(gen) && kt2>kt2min) {
            double maxlogmed = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                if (q2>kt2) continue;
                double phasemax = inte_C((2*zmax*E)/(t*kt2min)) - inte_C((2*zmin*E)/(t*kt2));
                maxlogmed += q2 * phasemax ;
            }
            maxlogmed *= 2. / kt2min * CAoverCR;
            double Crad = CR_2overb0 * (logvac + maxlogmed);
            double r = flat_gen(gen);
            kt2 = mu2 * std::pow(kt2/mu2, std::pow(r, 1.0/Crad) );
            double logmed = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                if (q2>kt2) continue;
                double phase = inte_C((2*zmax*E)/(t*kt2)) - inte_C((2*zmin*E)/(t*kt2));
                logmed += q2 * phase;
            }
            logmed *= 2./kt2*CAoverCR;
            acceptance = (logvac + logmed) / (logvac + maxlogmed);
        }
    } else {
        // Genearlized formula
        double maxdiffz = 1./zmin - 1./zmax + 2.*logvac;
        while (acceptance<flat_gen(gen) && kt2>kt2min) {
            double maxmedcoeff = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                maxmedcoeff += t*(kt2 + q2);
            }
            maxmedcoeff *= CAoverCR/(2.*E)*maxdiffz;
            double Crad = CR_2overb0 * (logvac + maxmedcoeff);
            double r = flat_gen(gen);
            kt2 = mu2 * std::pow(kt2/mu2, std::pow(r, 1.0/Crad) );
            // compute acceptance
            double medcoeff = 0.;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                medcoeff += t*(kt2 + q2);
            }
            medcoeff *= CAoverCR/(2.*E)*maxdiffz;
            acceptance = (logvac + medcoeff) / (logvac + maxmedcoeff);
        }
    }
    return (kt2>kt2min);
}

double eHIJING::sample_z_stochastic(double & z, int pid,
                        double E,
                        double kt2,
                        std::vector<double> qt2s,
                        std::vector<double> ts,
                        std::vector<double> phis) {
    double CR = (pid==21)? CA : CF;
    double zmin = std::min(.4/E, .4);
    double zmax = 1. - zmin;

    int Ncolls = ts.size();
    double acceptance = 0.;
    double weight = 1.0;
    if (mode_==0){
        while (acceptance<flat_gen(gen)) {
            z = zmin*std::pow(zmax/zmin, flat_gen(gen));
            double tauf = 2.*z*(1.0-z)*E/kt2;
            double w = 1.0, wmax = 1.0;
            for (int i=0; i<Ncolls; i++){
                double q2 = qt2s[i], t = ts[i];
                if (q2>kt2) continue;
                w += 2.*q2/kt2 * CA / CR * (1.-cos(t/tauf));
                wmax += 4.*q2/kt2 * CA / CR;
            }
            acceptance = w/wmax;
        }
    } else {
        double a = 0.;
        for (int i=0; i<Ncolls; i++){
            double q2 = qt2s[i], t = ts[i];
            a += (kt2 +q2)*t;
        }
        a *= CA/CR/(2.*E);
        double Norm = (1. + 2.*a)*std::log(zmax/zmin) + a*(1./zmin-1/zmax);
        // sample z ~ 1/z + a/[z^2(1-z)]
      	double left = zmin;
      	double right = zmax;
      	double mid = (left+right)/2.;
      	double r = flat_gen(gen);
      	while(right-left > 1e-3) {
        		mid = (left+right)/2.;
        		double fmid = ( (1+a)*std::log(mid/zmin)
        			    + a*(1./zmin-1/mid)
        			    + a*std::log(zmax/(1-mid)) ) / Norm;
        		if (fmid<r) left = mid;
        		else right = mid;
      	}
        z = mid;
        if (Ncolls==0) weight=1.;
        else {
           double wmax = CR/CA;
           for (int i=0; i<Ncolls; i++){
               double q2 = qt2s[i], t = ts[i];
               wmax += (kt2+q2)*t;
           }
           wmax /= (2.*z*(1-z)*E);

           double w = CR/CA;
           for (int i=0; i<Ncolls; i++){
               double dw;
               double q2 = qt2s[i], t = ts[i], phi = phis[i];
	       double A = (kt2+q2)*t/(2.*z*(1-z)*E),
                      B = 2.*std::sqrt(kt2*q2)*t/(2.*z*(1-z)*E);
               double X = B/A;
               double Y = A*(1.-X);
               if (X<.5){
                    double jv0 = std::cyl_bessel_j(0,B),
                           jv1 = std::cyl_bessel_j(1,B);
                    dw = - X*std::sin(A)*jv1
                         + X*X * ( .5 + X*std::cos(A) * (jv1/B - jv0) );
               } else {
                   if (Y>10.){
                       dw = 1./std::sqrt(1.-X*X)-1;
                   }
                   else {
                       dw = GHT_Angular_Table.interpolate({X, std::log(1.+Y)});
                   }
               }
               w += dw;
           }
           weight = w/wmax;
        }
    }
    return std::min(std::max(weight,0.),1.);
}

} //End eHIJING namespace