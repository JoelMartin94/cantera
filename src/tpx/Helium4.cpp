/**
 * @file Helium4.cpp representation of substance Helium4.
 *
 * Values and functions are from "Thermodynamic Properties in SI" by W.C.
 * Reynolds AUTHOR: me@rebeccahhunt.com: GCEP, Stanford University
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Helium4.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{

/*
 * Helium4 constants
 */


static const double Tmn = 2.177; // [K] minimum temperature for which calculations are valid
static const double Tmx = 1501.0; // [K] maximum temperature for which calculations are valid
static const double Tc=5.2014; // [K] critical temperature
static const double Roc=69.64; // [kg/m^3] critical density
static const double To=2.177; // [K] reference Temperature
static const double R=2077.22578699; // [] gas constant for CO2 J/kg/K
static const double Gamma1 = 1.56047072875E-4;
static const double Gamma2 = 3.12094145751E-5; // [??]
static const double Gamma3 = 3.12094145751E-5; // [??]
static const double u0=1.8712207E4; // [] internal energy at To
static const double s0=1.0812833E4; // [] entropy at To
static const double Tp=250; // [K] ??
static const double Pc=0.22746E6; // [Pa] critical pressure
static const double M=4.0026; // [kg/kmol] molar density

// array Ahel is used by the function named Pp
static const double Ahel1[]= {
	-2.63717841606E-4,		//1
	-5.79620044301E-2,		//2
	6.04727743809,			//3
	3.86500111589E1,		//4
	-2.75796664744E2, 		//5
	-4.96960774707E2,		//6
	2.04341052964E3,		//7
	-2.66595676810E3,		//8
	1.07968703317E3,		//9
	2.33740311250E-1,		//10
	-5.14034417722,			//11
	3.08419481342E1, 		//12
	-1.67047385071E2, 		//13
	5.24045883077E2, 		//14
	-8.07915654647E2,		//15
	6.31099960781E2,		//16
	-2.45791511511E2,		//17
	1.47668657398E-2,		//18
	-2.53062442742E-1,		//19
	7.33463898526E-1, 		//20
	2.92163822280E-1,		//21
	4.07953759561E-3,		//22
	-3.73905300971E-2,		//23
	1.36171997779E-1,		//24
	-2.47415495892E-1,		//25
	2.33727221372E-1,		//26
	-9.44142746383E-2		//27
	-3.72006192405E-6,		//28
	1.59283523218E-5,		//29
	7.75248537108,			//30
	-4.13169817472E1,		//31
	5.40743659299E1,		//32
	5.34172600153E-4,		//33
	-1.05413018834E-3,		//34
	-8.82580260817E-4		//35
};

static const double Ahel2[]= {
	-2.63717841606E-4,		//1
	-5.79620044301E-2,		//2
	6.04727743809,			//3
	3.86500111589E1,		//4
	-2.75796664744E2, 		//5
	-4.96960774707E2,		//6
	2.04341052964E3,		//7
	-2.66595676810E3,		//8
	1.07968703317E3,		//9
	3.23316248529E-2, 		//10
	2.01417823467,			//11
	-3.20336592218E1,		//12
	1.17952847254E2,		//13
	-2.72064513304E2,		//14
	8.06705554799E2,		//15
	-6.34863771449E2,		//16
	4.23944026969E2,		//17
	-1.26804959063E-2,		//18
	5.58960362485E-2,		//19
	5.81328684698E-1,		//20
	-1.03365680210,			//21
	-1.01057001312E-3,		//22
	8.40859671873E-3,		//23
	-2.48181422872E-2,		//24
	3.24270326025E-2,		//25
	-1.04566294786E-2,		//26
	-1.05412341221E-2,		//27
	-1.04201749588E-6,		//28
	1.09726080203E-5,		//29
	1.24933778088E1,		//30
	-1.41252424541E2,		//31
	-2.38228039845E2,		//32
	2.65139980533E-4,		//33
	3.33310756017E-3,		//34
	-2.41601688592E-3		//35
};

static const double Ahel3[]= {
	-2.63717841606E-4,		//1
	-5.79620044301E-2,		//2
	6.04727743809,			//3
	3.86500111589E1,		//4
	-2.75796664744E2, 		//5
	-4.96960774707E2,		//6
	2.04341052964E3,		//7
	-2.66595676810E3,		//8
	1.07968703317E3,		//9
	-5.69281410539E-2,		//10
	2.54082433493,			//11
	-4.33612764494E1,		//12
	2.32901880818E2,		//13
	-6.88289870860E2,		//14
	2.12493828516E3,		//15
	-2.69258356337E3,		//16
	1.42625846393E3,		//17
	7.76178949940E-4,		//18
	6.75967782095E-2,		//19
	9.09992115812E-2,		//20
	-3.81211874106E-1,		//21
	-2.30068006523E-5,		//22
	4.02950826349E-5,		//23
	1.07511466109E-3,		//24
	-4.93747339170E-3,		//25
	1.11576934297E-2,		//26
	-1.23679512941E-2,		//27
	-3.64745210287E-7,		//28
	1.02807881652E-5,		//29
	8.98703364016,			//30
	-2.28140026278E2,		//31
	5.33588707469,			//32
	1.06067862115E-4,		//33
	-4.46441499497E-3,		//34
	3.80683087199E-3		//35
};

const double *Ahel = Ahel1;
double Gamma = Gamma1;

// array F is used by the function named Psat
static const double F[]= {
	-3.9394635287, 		//1
	1.3925998798E2,		//2
	-1.6407741565E3,	//3
	1.1974557102E4,		//4
	-5.5283309818E4,	//5
	1.16621956504E5,	//6
	-3.2521282840E5,	//7
	3.9884322750E5,		//8
	-2.771806992E5,		//9
	8.3395204183E4		//10
};

// array D is used by the function ldens
static const double D[]= {
    6.6940000000E1,		//1
	1.2874326484E2,		//2
	-4.3128217346E2,	//3
	1.7851911824E3,		//4
	-3.3509624489E3,	//5
	3.0344215824E3,		//6
	-1.0981289602E3		//7
};

// double G is used by the function sp
static const double G = 3115.85;

double Helium4::C(int j,double Tinverse, double T2inverse, const double *Ahel)
{
	double sum=0.0;
	int i;

    switch (j) {
    case 0:
    	for (i=0; i<=8; i++) {
    	        sum += *(Ahel + i)*pow(T, 2.-0.5*double(i));
    	    }
        return sum;
    case 1:
    	for (i=9; i<=16; i++) {
    	        sum += *(Ahel + i)*pow(T, 5.5-0.5*double(i));
    	    }
        return sum;
    case 2:
    	for (i=17; i<=20; i++) {
    	        sum += *(Ahel + i)*pow(T, 17.5-double(i));
    	    }
        return sum;
    case 3:
    	for (i=21; i<=26; i++) {
    	        sum += *(Ahel + i)*pow(T, 5.75-0.25*double(i));
    	    }
        return sum;
    case 4:
        return *(Ahel + 27) +
        	   *(Ahel + 28)*Tinverse;
    case 5:
        return *(Ahel + 29) +
        		*(Ahel + 30) *Tinverse +
				*(Ahel + 31) *T2inverse;
    case 6:
        return *(Ahel + 32) +
        		*(Ahel + 33) *Tinverse +
        				*(Ahel + 34) *T2inverse;
    default:
    	throw CanteraError("Helium4::C",
    	    	                           "Index out of range. j = {}", j);
    }
}

double Helium4::Cprime(int j, double T2inverse, double T3inverse, double T4inverse, const double *Ahel)
{
	double sum = 0.0;
	int i;

    switch (j) {
    case 0:
    	for (i=0; i<=8; i++) {
    	        sum += *(Ahel + i)*(2.-0.5*double(i))*pow(T, 1.-0.5*double(i));
    	    }
        return sum;
    case 1:
    	for (i=9; i<=16; i++) {
    	        sum += *(Ahel + i)*(5.5-0.5*double(i))*pow(T, 4.5-0.5*double(i));
    	    }
        return sum;
    case 2:
    	for (i=17; i<=20; i++) {
    	        sum += *(Ahel + i)*(17.5-double(i))*pow(T, 16.5-double(i));
    	    }
        return sum;
    case 3:
    	for (i=21; i<=26; i++) {
    	        sum += *(Ahel + i)*(5.75-0.25*double(i))*pow(T, 4.75-0.25*double(i));
    	    }
        return sum;
    case 4:
        return -*(Ahel + 28)*T2inverse;
    case 5:
        return
            - *(Ahel + 30) *T2inverse +
            - 2.0* *(Ahel + 31) *T3inverse;
    case 6:
        return
            - *(Ahel + 33) *T2inverse +
            - 2.0* *(Ahel + 34) *T3inverse;
    default:
    	throw CanteraError("Helium4::Cprime",
    	    	                           "Index out of range. j = {}", j);
    }
}

double Helium4::I(int j, double ergho, double Gamma)
{
    switch (j) {
    case 0:
        return Rho;
    case 1:
        return pow(Rho, 2)/2;
    case 2:
        return pow(Rho, 3)/ 3;
    case 3:
        return pow(Rho, 4)/ 4;
    case 4:
        return pow(Rho, 5)/ 5;
    case 5:
        return (1 - ergho) / double(2 * Gamma);
    case 6:
        return (1 - ergho * double(Gamma * pow(Rho,2) + double(1)))/ double(2 * Gamma * Gamma);
    default:
    	throw CanteraError("Helium4::I",
    	                           "Index out of range. j = {}", j);
    	//return 0.0;
    }
}

double Helium4::H(int i, double egrho)
{
    if (i <= 4) {
        return pow(Rho,i+2);
    } else if (i == 5) {
        return pow(Rho,3)*egrho;
    } else if (i == 6) {
        return pow(Rho,5)*egrho;
    } else {
        return 0;
    }
}

double Helium4::up()
{
    double Tinverse = 1.0/T;
    double T2inverse = pow(T, -2);
    double T3inverse = pow(T, -3);
    double T4inverse = pow(T, -4);
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    // Equation C-6 integrated
    sum += G*(T-To);
    int i;
    for (i=0; i<=6; i++) {
        sum += I(i,egrho, Gamma) *
               (C(i, Tinverse, T2inverse, Ahel) - T*Cprime(i,T2inverse, T3inverse, T4inverse, Ahel));
    }
    sum += u0;
    return sum + m_energy_offset;
}

double Helium4::sp()
{
    double T2inverse = pow(T, -2);
    double T3inverse = pow(T, -3);
    double T4inverse = pow(T, -4);
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    sum += G*log(T/To);
    for (int i=0; i<=6; i++) {
        sum -= Cprime(i,T2inverse, T3inverse, T4inverse, Ahel)*I(i,egrho,Gamma);
    }
    sum += s0 - R*log(Rho);
    return sum + m_entropy_offset;
}

double Helium4::ldens()
{
    double xx=1-(T/Tc), sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Helium4::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=6; i++) {
        sum+=D[i-1]*pow(xx,double(i-1)/3.0);
    }
    return sum;
}

double Helium4::Pp()
{
    double Tinverse = pow(T,-1);
    double T2inverse = pow(T, -2);
    double P = Rho*R*T;

    if (2.0 <= T && T<=Tc)
    {
    	if (Rho <= ldens())
    	{
        	//Region I
    		Ahel = Ahel1;
    		Gamma = Gamma1;
    	}
    	else
    	{
    		//Region II
    		Ahel = Ahel2;
    		Gamma = Gamma2;
    	}
    }
    else if (Tc <= T && T <= 10.0)
    {
    	if (Rho <= Roc)
    	{
    		//Region I
    		Ahel = Ahel1;
    		Gamma = Gamma1;
    	}
    	else
    	{
        	//Region II
    		Ahel = Ahel2;
    		Gamma = Gamma2;
    	}
    }
    else if (10 < T && T < 15.0)
    {
    	throw CanteraError("Helium4::Pp",
    	    	    	                           "Region not implemented. T = {}", T, " Rho = {}", Rho);
    	if (Rho <= Roc)
    	{
    		//Region I + Region III
    	}
    	else
    	{
        	//Region II + Region III
    	}
    }
    else if (15.0 <= T)
    {
    	//Region III
    	Ahel = Ahel3;
    	Gamma = Gamma3;
    }
    else
    {
    	throw CanteraError("Helium4::Pp",
    	    	                           "Invalid region for. T = {}", T, " Rho = {}", Rho);
    }

    double egrho = exp(-Gamma*Rho*Rho);


    // when i=0 we are on second sum of equation (where rho^2)
    for (int i=0; i<=6; i++) {
        P += C(i,Tinverse, T2inverse, Ahel)*H(i,egrho);
    }
    return P;
}

double Helium4::Psat()
{
    double sum=0,P;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Helium4::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=10; i++) {
        sum += F[i-1] * pow(T,double(2-i));
    }

    P=exp(sum);
    return P;
}


// The following functions allow users to get the properties of Helium4
// that are not dependent on the state

double Helium4::Tcrit()
{
    return Tc;
}
double Helium4::Pcrit()
{
    return Pc;
}
double Helium4::Vcrit()
{
    return 1.0/Roc;
}
double Helium4::Tmin()
{
    return Tmn;
}
double Helium4::Tmax()
{
    return Tmx;
}
double Helium4::MolWt()
{
    return M;
}

}
