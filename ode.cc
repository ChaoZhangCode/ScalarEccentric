#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_ellint.h>
#include <algorithm>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <complex>
#include <cmath>

#include "Interpolant.h"
#include "global.h"
#include "Utility.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <iomanip>      // std::setprecision
#include <cstring>

#include "dIdt8H_5PNe10.h"
#include "ode.hh"

#define pn5_Y
#define pn5_citation1 Pn5_citation
__deriv__
void pn5_base_func(double* pdot, double* edot, double* Ydot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double Y, double* additional_args)
{
    // evaluate ODEs

    // the frequency variables are pointers!
    double x = Y_to_xI(a, p, e, Y);
    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, x);

	int Nv = 10;
    int ne = 10;
    *pdot = epsilon * dpdt8H_5PNe10 (a, p, e, Y, Nv, ne);

    // needs adjustment for validity
    Nv = 10;
    ne = 8;
	*edot = epsilon * dedt8H_5PNe10 (a, p, e, Y, Nv, ne);

    Nv = 7;
    ne = 10;
    *Ydot = epsilon * dYdt8H_5PNe10 (a, p, e, Y, Nv, ne);

}


// Initialize flux data for inspiral calculations
void load_and_interpolate_flux_data(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "few/files/FluxNewMinusPNScaled_fixed_y_order.dat";
    fp = few_dir + fp;
	ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxNewMinusPNScaled_fixed_y_order.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ys, es, Edots, Ldots;
	double y, e, Edot, Ldot;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

		Flux_ss >> y >> e >> Edot >> Ldot;

		ys.push_back(y);
		es.push_back(e);
		Edots.push_back(Edot);
		Ldots.push_back(Ldot);
	}

	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ys.begin(), ys.end() );
	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

	sort( es.begin(), es.end() );
	es.erase( unique( es.begin(), es.end() ), es.end() );

	Interpolant *Edot_interp = new Interpolant(ys, es, Edots);
	Interpolant *Ldot_interp = new Interpolant(ys, es, Ldots);

	interps->Edot = Edot_interp;
	interps->Ldot = Ldot_interp;

}


// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFlux::SchwarzEccFlux(std::string few_dir)
{
    interps = new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place
    load_and_interpolate_flux_data(interps, few_dir);
	//load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}

#define SchwarzEccFlux_num_add_args 0
#define SchwarzEccFlux_spinless
#define SchwarzEccFlux_equatorial
#define SchwarzEccFlux_file1 FluxNewMinusPNScaled_fixed_y_order.dat
__deriv__
void SchwarzEccFlux::deriv_func(double* pdot, double* edot, double* xdot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double x, double* additional_args)
{
    if ((6.0 + 2. * e) > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    SchwarzschildGeoCoordinateFrequencies(Omega_phi, Omega_r, p, e);
    *Omega_theta = *Omega_phi;

    double y1 = log((p -2.*e - 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

	double yPN = pow((*Omega_phi),2./3.);

	double EdotPN = (96 + 292*Power(e,2) + 37*Power(e,4))/(15.*Power(1 - Power(e,2),3.5)) * pow(yPN, 5);
	double LdotPN = (4*(8 + 7*Power(e,2)))/(5.*Power(-1 + Power(e,2),2)) * pow(yPN, 7./2.);

	double Edot = -epsilon*(interps->Edot->eval(y1, e)*pow(yPN,6.) + EdotPN);

	double Ldot = -epsilon*(interps->Ldot->eval(y1, e)*pow(yPN,9./2.) + LdotPN);

	*pdot = (-2*(Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*(3 + Power(e,2) - p)*Power(p,1.5) + Ldot*Power(-4 + p,2)*Sqrt(-3 - Power(e,2) + p)))/(4*Power(e,2) - Power(-6 + p,2));

    // handle e = 0.0
	if (e > 0.)
    {
        *edot = -((Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*Power(p,1.5)*
            	  (18 + 2*Power(e,4) - 3*Power(e,2)*(-4 + p) - 9*p + Power(p,2)) +
            	 (-1 + Power(e,2))*Ldot*Sqrt(-3 - Power(e,2) + p)*(12 + 4*Power(e,2) - 8*p + Power(p,2)))/
            	(e*(4*Power(e,2) - Power(-6 + p,2))*p));
    }
    else
    {
        *edot = 0.0;
    }

    *xdot = 0.0;
}


// destructor
SchwarzEccFlux::~SchwarzEccFlux()
{

    delete interps->Edot;
    delete interps->Ldot;
    delete interps;


}













void load_and_interpolate_flux_data2(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "few/files/FluxNewMinusPNScaled_fixed_y_order.dat";
    fp = few_dir + fp;
	ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxNewMinusPNScaled_fixed_y_order.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ys, es, Edots, Ldots;
	double y, e, Edot, Ldot;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

		Flux_ss >> y >> e >> Edot >> Ldot;

		ys.push_back(y);
		es.push_back(e);
		Edots.push_back(Edot);
		Ldots.push_back(Ldot);
	}

	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ys.begin(), ys.end() );
	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

	sort( es.begin(), es.end() );
	es.erase( unique( es.begin(), es.end() ), es.end() );

	Interpolant *Edot_interp = new Interpolant(ys, es, Edots);
	Interpolant *Ldot_interp = new Interpolant(ys, es, Ldots);

	interps->Edot = Edot_interp;
	interps->Ldot = Ldot_interp;

}

void load_and_interpolate_scalar_eccentric_flux_data(struct interp_params *interps, const std::string& few_dir){

	// Load and interpolate the flux data
    std::string fp = "few/files/FluxScalarpefinal.txt";
    fp = few_dir + fp;
	ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxNewMinusPNScaled_fixed_y_order.dat did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

	// Load the flux data into arrays
	string Flux_string;
	vector<double> ys, es, Edots, Ldots;
	double y, e, Edot, Ldot;
	while(getline(Flux_file, Flux_string)){

		stringstream Flux_ss(Flux_string);

		Flux_ss >> y >> e >> Edot >> Ldot;

		ys.push_back(y);
		es.push_back(e);
		Edots.push_back(Edot);
		Ldots.push_back(Ldot);
	}

	// Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)
	sort( ys.begin(), ys.end() );
	ys.erase( unique( ys.begin(), ys.end() ), ys.end() );

	sort( es.begin(), es.end() );
	es.erase( unique( es.begin(), es.end() ), es.end() );

	Interpolant *Edot_interp = new Interpolant(ys, es, Edots);
	Interpolant *Ldot_interp = new Interpolant(ys, es, Ldots);

	interps->Edot = Edot_interp;
	interps->Ldot = Ldot_interp;

}


// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFluxwithscalar::SchwarzEccFluxwithscalar(std::string few_dir)
{
    interps = new interp_params;
    interpsscalar= new interp_params;

    load_and_interpolate_flux_data2(interps, few_dir);
    load_and_interpolate_scalar_eccentric_flux_data(interpsscalar, few_dir);

}

#define SchwarzEccFluxwithscalar_num_add_args 1
#define SchwarzEccFluxwithscalar_spinless
#define SchwarzEccFluxwithscalar_equatorial
#define SchwarzEccFluxwithscalar_file1 FluxNewMinusPNScaled_fixed_y_order.dat
#define SchwarzEccFluxwithscalar_file2 FluxScalarpefinal.txt
__deriv__
void SchwarzEccFluxwithscalar::deriv_func(double* pdot, double* edot, double* xdot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double x, double* additional_args)
{
    if ((6.0 + 2. * e) > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }

    double scharge=additional_args[0];

    SchwarzschildGeoCoordinateFrequencies(Omega_phi, Omega_r, p, e);
    *Omega_theta = *Omega_phi;

    double y1 = log((p -2.*e - 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

	double yPN = pow((*Omega_phi),2./3.);

	double EdotPNtensor = (96 + 292*Power(e,2) + 37*Power(e,4))/(15.*Power(1 - Power(e,2),3.5)) * pow(yPN, 5);
	double LdotPNtensor  = (4*(8 + 7*Power(e,2)))/(5.*Power(-1 + Power(e,2),2)) * pow(yPN, 7./2.);

	double Edottensor = -epsilon*(interps->Edot->eval(y1, e)*pow(yPN,6.) + EdotPNtensor);

	double Ldottensor = -epsilon*(interps->Ldot->eval(y1, e)*pow(yPN,9./2.) + LdotPNtensor);

	double EdotPNscalar = 1./3.* pow(yPN, 4)/pow((1.-Power(e,2)),2.5)*(1+0.5*Power(e,2));

	double LdotPNscalar  =1./3.* pow(yPN, 2.5)/(1.-Power(e,2));

        double Edotscalar = -epsilon*(interpsscalar->Edot->eval(y1,e)*pow(yPN,5.)+ EdotPNscalar);

        double Ldotscalar=-epsilon*(interpsscalar->Ldot->eval(y1,e)*pow(yPN,3.5)+ LdotPNscalar);



        double Edot =Edottensor+Edotscalar*pow(scharge,2.);

        double Ldot=Ldottensor+Ldotscalar*pow(scharge,2.);

	*pdot = (-2*(Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*(3 + Power(e,2) - p)*Power(p,1.5) + Ldot*Power(-4 + p,2)*Sqrt(-3 - Power(e,2) + p)))/(4*Power(e,2) - Power(-6 + p,2));

    // handle e = 0.0
	if (e > 0.)
    {
        *edot = -((Edot*Sqrt((4*Power(e,2) - Power(-2 + p,2))/(3 + Power(e,2) - p))*Power(p,1.5)*
            	  (18 + 2*Power(e,4) - 3*Power(e,2)*(-4 + p) - 9*p + Power(p,2)) +
            	 (-1 + Power(e,2))*Ldot*Sqrt(-3 - Power(e,2) + p)*(12 + 4*Power(e,2) - 8*p + Power(p,2)))/
            	(e*(4*Power(e,2) - Power(-6 + p,2))*p));
    }
    else
    {
        *edot = 0.0;
    }

    *xdot = 0.0;
}


// destructor
SchwarzEccFluxwithscalar::~SchwarzEccFluxwithscalar()
{

    delete interps->Edot;
    delete interps->Ldot;
    delete interps;

    delete interpsscalar->Edot;
    delete interpsscalar->Ldot;
    delete interpsscalar;

}





































// Initialize flux data for inspiral calculations
void load_and_interpolate_tensor_flux_data(struct interp_params *interps, const std::string& few_dir){

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxTensor.txt";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxTensor.txt did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys,  Edots;
    double y,  Edot;
    while(getline(Flux_file, Flux_string)){

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y  >> Edot ;

        ys.push_back(y);
        Edots.push_back(Edot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)


    Interpolant *Edot_interp = new Interpolant(ys, Edots);

    interps->Edot = Edot_interp;

}

// Initialize flux data for inspiral calculations
void load_and_interpolate_vector_flux_data(struct interp_params *interps, const std::string& few_dir){

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxVector.txt";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxVector.txt did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys,  Edots;
    double y,  Edot;
    while(getline(Flux_file, Flux_string)){

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y  >> Edot ;

        ys.push_back(y);
        Edots.push_back(Edot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)


    Interpolant *Edot_interp = new Interpolant(ys, Edots);

    interps->Edot = Edot_interp;

}

// Initialize flux data for inspiral calculations
void load_and_interpolate_scalar_flux_data(struct interp_params *interps, const std::string& few_dir){

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxScalar.txt";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxScalar.txt did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys,  Edots;
    double y,  Edot;
    while(getline(Flux_file, Flux_string)){

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y  >> Edot ;

        ys.push_back(y);
        Edots.push_back(Edot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)


    Interpolant *Edot_interp = new Interpolant(ys, Edots);

    interps->Edot = Edot_interp;

}


// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFlux2::SchwarzEccFlux2(std::string few_dir)
{
    interpstensor = new interp_params;
    interpsvector = new interp_params;
    interpsscalar= new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place

    load_and_interpolate_tensor_flux_data(interpstensor, few_dir);
    load_and_interpolate_vector_flux_data(interpsvector, few_dir);
    load_and_interpolate_scalar_flux_data(interpsscalar, few_dir);
    //load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}

#define SchwarzEccFlux2_num_add_args 2
#define SchwarzEccFlux2_spinless
#define SchwarzEccFlux2_equatorial
#define SchwarzEccFlux2_file1 FluxScalar.txt
#define SchwarzEccFlux2_file2 FluxTensor.txt
#define SchwarzEccFlux2_file3 FluxVector.txt
__deriv__
void SchwarzEccFlux2::deriv_func(double* pdot, double* edot, double* xdot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double x, double* additional_args)
{
    if (6.0  > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *xdot = 0.0;
        return;
    }
    double scharge=additional_args[0];
    double vcharge=additional_args[1];

    SchwarzschildGeoCoordinateFrequencies(Omega_phi, Omega_r, p, e);
    *Omega_theta = *Omega_phi;

    double y1 = log((p- 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

	double yPN = pow((*Omega_phi),2./3.);

	double EdotPNtensor = 32./5. * pow(yPN, 5);

	double EdotPNvector = 1./3. * pow(yPN, 4);

	double Edottensor = -epsilon*(interpstensor->Edot->eval(y1)*pow(yPN,6.) + EdotPNtensor);


	double Edotvector = -epsilon*(interpsvector->Edot->eval(y1)*pow(yPN,5.)+ EdotPNvector);


        double Edotscalar = -epsilon*(interpsscalar->Edot->eval(y1)*pow(yPN,5.)+ EdotPNvector);



    double Edot =Edottensor+Edotvector*Power(vcharge,2)+Edotscalar*Power(scharge,2);


    *pdot = 1.0/((p-6.0)/2.0/Power(p*(p-3.0),1.5))*Edot;
    *edot = 0.0;
    *xdot = 0.0;
    
}


// destructor
SchwarzEccFlux2::~SchwarzEccFlux2()
{
   delete interpsscalar->Edot;
    delete interpsscalar;

    delete interpstensor->Edot;
    delete interpstensor;

    delete interpsvector->Edot;
    delete interpsvector;

}





// Initialize flux data for inspiral calculations
void load_and_interpolate_tensor_flux_datakerr(struct interp_params *interps, const std::string& few_dir){

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxTensorkerr09.txt";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxTensorkerr.txt did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys,  Edots;
    double y,  Edot;
    while(getline(Flux_file, Flux_string)){

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y  >> Edot ;

        ys.push_back(y);
        Edots.push_back(Edot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)


    Interpolant *Edot_interp = new Interpolant(ys, Edots);

    interps->Edot = Edot_interp;

}

// Initialize flux data for inspiral calculations
void load_and_interpolate_scalar_flux_datakerr(struct interp_params *interps, const std::string& few_dir){

    // Load and interpolate the flux data
    std::string fp = "few/files/FluxScalarkerr09.txt";
    fp = few_dir + fp;
    ifstream Flux_file(fp);

    if (Flux_file.fail())
    {
        throw std::runtime_error("The file FluxScalarkerr09.txt did not open sucessfully. Make sure it is located in the proper directory (Path/to/Installation/few/files/).");
    }

    // Load the flux data into arrays
    string Flux_string;
    vector<double> ys,  Edots;
    double y,  Edot;
    while(getline(Flux_file, Flux_string)){

        stringstream Flux_ss(Flux_string);

        Flux_ss >> y  >> Edot ;

        ys.push_back(y);
        Edots.push_back(Edot);
    }

    // Remove duplicate elements (only works if ys are perfectly repeating with no round off errors)


    Interpolant *Edot_interp = new Interpolant(ys, Edots);

    interps->Edot = Edot_interp;

}

// Class to carry gsl interpolants for the inspiral data
// also executes inspiral calculations
SchwarzEccFlux3::SchwarzEccFlux3(std::string few_dir)
{
    interpstensor = new interp_params;
    interpsscalar = new interp_params;

    // prepare the data
    // python will download the data if
    // the user does not have it in the correct place
    load_and_interpolate_tensor_flux_datakerr(interpstensor, few_dir);
    load_and_interpolate_scalar_flux_datakerr(interpsscalar, few_dir);
    //load_and_interpolate_amp_vec_norm_data(&amp_vec_norm_interp, few_dir);
}
#define SchwarzEccFlux3_num_add_args 1
#define SchwarzEccFlux3_Y
#define SchwarzEccFlux3_file1 FluxScalarkerr09.txt
#define SchwarzEccFlux3_file2 FluxTensorkerr09.txt
__deriv__
void SchwarzEccFlux3::deriv_func(double* pdot, double* edot, double* Ydot,
                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                  double epsilon, double a, double p, double e, double Y, double* additional_args)
{

    if (get_separatrix(a,e,1.0) > p)
    {
        *pdot = 0.0;
        *edot = 0.0;
        *Ydot = 0.0;
        return;
    }


    KerrGeoCoordinateFrequencies(Omega_phi, Omega_theta, Omega_r, a, p, e, 1.0);


    double scharge=additional_args[0];


    double y1 = log((p- 2.1));

    // evaluate ODEs, starting with PN contribution, then interpolating over remaining flux contribution

	double yPN = pow((*Omega_phi),2./3.);

	double EdotPNtensor = 32./5. * pow(yPN, 5)+8./15.*pow(yPN,6.5)*a*(-49);

	double EdotPNscalar = 1./3. * pow(yPN, 4)*(1-2*yPN+(2*Pi-4*a)*pow(yPN,1.5));

	double Edottensor = -epsilon*(interpstensor->Edot->eval(y1)*pow(yPN,6.) + EdotPNtensor);


	double Edotscalar = -epsilon*(interpsscalar->Edot->eval(y1)*pow(yPN,5.)+ EdotPNscalar);



        double Edot =Edottensor+Edotscalar*Power(scharge,2);



    *pdot=1./(-3*Power(a,2)+8*a*Power(p,0.5)+(p-6)*p)*(2*Power(2*a+(-3+p)*Power(p,0.5),1.5)*Power(p,7.0/4.0))*Edot;

    *edot = 0.0;
    *Ydot = 0.0;
    
}


// destructor
SchwarzEccFlux3::~SchwarzEccFlux3()
{

    delete interpstensor->Edot;
    delete interpstensor;

    delete interpsscalar->Edot;
    delete interpsscalar;

}

                pn5::pn5(std::string few_dir){}

                pn5::~pn5(){}

                void pn5::deriv_func(double* pdot, double* edot, double* Ydot,
                                  double* Omega_phi, double* Omega_theta, double* Omega_r,
                                  double epsilon, double a, double p, double e, double Y, double* additional_args)
                {
                    pn5_base_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                  epsilon, a, p, e, Y, additional_args);
                }
            

    ODECarrier::ODECarrier(std::string func_name_, std::string few_dir_)
    {
        func_name = func_name_;
        few_dir = few_dir_;
    
            if (func_name == "pn5")
            {
        
                pn5* temp = new pn5(few_dir);

                func = (void*) temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux")
            {
        
                SchwarzEccFlux* temp = new SchwarzEccFlux(few_dir);

                func = (void*) temp;

            
            }

        
            else if (func_name == "SchwarzEccFluxwithscalar")
            {
        
                SchwarzEccFluxwithscalar* temp = new SchwarzEccFluxwithscalar(few_dir);

                func = (void*) temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux2")
            {
        
                SchwarzEccFlux2* temp = new SchwarzEccFlux2(few_dir);

                func = (void*) temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux3")
            {
        
                SchwarzEccFlux3* temp = new SchwarzEccFlux3(few_dir);

                func = (void*) temp;

            
            }

        
    }
    

    void ODECarrier::get_derivatives(double* pdot, double* edot, double* Ydot,
                      double* Omega_phi, double* Omega_theta, double* Omega_r,
                      double epsilon, double a, double p, double e, double Y, double* additional_args)
    {
    
            if (func_name == "pn5")
            {
        
                pn5* temp = (pn5*)func;

                temp->deriv_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                epsilon, a, p, e, Y, additional_args);

            
            }

        
            else if (func_name == "SchwarzEccFlux")
            {
        
                SchwarzEccFlux* temp = (SchwarzEccFlux*)func;

                temp->deriv_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                epsilon, a, p, e, Y, additional_args);

            
            }

        
            else if (func_name == "SchwarzEccFluxwithscalar")
            {
        
                SchwarzEccFluxwithscalar* temp = (SchwarzEccFluxwithscalar*)func;

                temp->deriv_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                epsilon, a, p, e, Y, additional_args);

            
            }

        
            else if (func_name == "SchwarzEccFlux2")
            {
        
                SchwarzEccFlux2* temp = (SchwarzEccFlux2*)func;

                temp->deriv_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                epsilon, a, p, e, Y, additional_args);

            
            }

        
            else if (func_name == "SchwarzEccFlux3")
            {
        
                SchwarzEccFlux3* temp = (SchwarzEccFlux3*)func;

                temp->deriv_func(pdot, edot, Ydot, Omega_phi, Omega_theta, Omega_r,
                                epsilon, a, p, e, Y, additional_args);

            
            }

        
    }
    

    ODECarrier::~ODECarrier()
    {
    
            if (func_name == "pn5")
            {
        
                pn5* temp = (pn5*)func;

                delete temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux")
            {
        
                SchwarzEccFlux* temp = (SchwarzEccFlux*)func;

                delete temp;

            
            }

        
            else if (func_name == "SchwarzEccFluxwithscalar")
            {
        
                SchwarzEccFluxwithscalar* temp = (SchwarzEccFluxwithscalar*)func;

                delete temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux2")
            {
        
                SchwarzEccFlux2* temp = (SchwarzEccFlux2*)func;

                delete temp;

            
            }

        
            else if (func_name == "SchwarzEccFlux3")
            {
        
                SchwarzEccFlux3* temp = (SchwarzEccFlux3*)func;

                delete temp;

            
            }

        
    }
    