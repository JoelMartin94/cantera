//! @file Helium4.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_HELIUM4_H
#define TPX_HELIUM4_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of helium4. Values and functions are
//! from "Thermodynamic Properties in SI" by W.C. Reynolds
class Helium4 : public Substance
{
public:
    Helium4() {
        m_name="helium4";
        m_formula="He";
    }

    double MolWt();
    double Tcrit();
    double Pcrit();
    double Vcrit();
    double Tmin();
    double Tmax();

    /*! Pressure. Equation P-5 from table 2.1 on page 122 in Reynolds. P(rho, T).
    	\f{eqnarray*}
	 	  P &=& \rho RT + \sum_{i=0}^{N-1} C_i(T)H_i(\rho)\\
	 	  	&=& \rho R T + \rho^2\sum_{i=0}^8 A_iT^{\frac{4-i}{2}} \\
	 	    && 	+ \rho^3 \sum_{i=9}^{16} A_iT^{\frac{11-i}{2}}\\
	 	    &&	+ \rho^4 \sum_{i=17}^{20} A_iT^{\frac{35-2i}{2}}\\
	 	    &&	+ \rho^5 \sum_{i=21}^{26} A_iT^{\frac{23-i}{4}}\\
	 	    &&	+ \rho^6 (A_{27} + A_{28}/T) \\
	 	    && 	+
	 	  	  \left[
	 	  	  	  \rho^3(A_{29} + A_{30}/T + A_{31}/T^2) +
	 	  	  	  \rho^5(A_{32} + A_{33}/T + A_{34}/T^2)
	 	  	  \right]e^{-\gamma\rho^2}
	 	\f}
	* Gamma exponent unsure tow or three
    */
    double Pp();

    /*! internal energy. See Reynolds eqn (15) on page 120 in section 2
      	\f[
	 	  u = \int_{T_0}^T c_v^0(T)dT + \sum_{i=0}^{N-1} [C_i -TC_i']I_i + u_0
	 	\f]
     */
    double up();

    /*! entropy. See Reynolds eqn (16) on page 120 in section 2
      	\f[
	 	  s = \int_{T_0}^T \frac{c_v^0(T)}{T}dT - R\log \rho - \sum_{i=0}^{N-1} C_i'I_i + s_0
	 	\f]
     */
    double sp();

    /*! Pressure at Saturation. Equation S-5 from table 2.3 on page 125 in Reynolds.
		\f[
		  \log P = \sum_{i=0}^{9} F_iT^{3-i}
		\f]
    */
    double Psat();

private:
    //! Liquid density. Equation D2 in Reynolds.
    double ldens();

    /*!
     * C returns a multiplier in each term of the sum in P-5, used in
     * conjunction with C in the function Pp
     * - j is used to represent which of the values in the summation to calculate
     * - j=0 is the second additive in the formula in reynolds
     * - j=1 is the third...
     * (this part does not include the multiplier rho^n)
     *
         \f{eqnarray*}
	 	  C_{i=0}^8 		&=& A_iT^{\frac{4-i}{2}}\\
	 	  C_{i=9}^{16}  	&=& A_iT^{\frac{11-i}{2}}\\
	 	  C_{i=17}^{20}  	&=&	A_iT^{\frac{35-2i}{2}}\\
	 	  C_{i=21}^{26}   	&=&	A_iT^{\frac{23-i}{4}}\\
	 	  C_{27}  			&=&	A_{27} \\
	 	  C_{28}  			&=& \frac{A_{28}}{T}\\
	 	  C_{29}  			&=& A_{29}\\
	 	  C_{30}  			&=& \frac{A_{30}}{T}\\
	 	  C_{31}  			&=& \frac{A_{31}}{T^2}\\
	 	  C_{32}  			&=& A_{32}\\
	 	  C_{33}  			&=& \frac{A_{33}}{T}\\
	 	  C_{34}  			&=& \frac{A_{34}}{T^2}
	 	\f}
     */
    double C(int jm, double, double);

    /*! Derivative of C(i)
     	 \f{eqnarray*}{
     	 	 C_i'				& = &	\frac{\partial C_i}{\partial T}\\
     	 	 {C'}_{i=0}^{8}		& = & A_{i} \left(\frac{4-i}{2}\right) T^{\frac{2-i}{2}} \\
     	 	 {C'}_{i=9}^{16} 	& = & A_{i} \left(\frac{11-i}{2}\right) T^{\frac{9-i}{2}}\\
     	 	 {C'}_{i=17}^{20} 	& = & A_{i} \left(\frac{35-2i}{2}\right) T^{\frac{33-2i}{2}}\\
     	 	 {C'}_{i=21}^{26} 	& = & A_{i} \left(\frac{23-i}{4}\right) T^{\frac{19-i}{4}}\\
     	 	 C'_{27}		& = & 0\\
     	 	 C'_{28} 	 	& = & -\frac{A_{28}}{T^{2}} \\
     	 	 C'_{29}		& = & 0\\
     	 	 C'_{30}		& = & -\frac{A_{30}}{T^{2}} \\
     	 	 C'_{31}		& = & -\frac{2 A_{31}}{T^{3}}\\
     	 	 C'_{32}		& = & 0\\
     	 	 C'_{33}		& = & - \frac{A_{33}}{T^{2}}\\
     	 	 C'_{34}		& = & - \frac{2 A_{34}}{T^{3}}
     	 \f}
     */
    double Cprime(int i, double, double, double);

    /*!
     * I = integral from o-rho { 1/(rho^2) * H(i, rho) d rho }
     * ( see eqaution 14b on page 120 in section 2 of Reynolds TPSI )
    	\f{eqnarray*}{
    		I_i 			&=& \int_0^{\rho} \frac{1}{\rho^2}H_i d \rho\\
    		I_{i=0}^8 		&=& \rho \\
		  	I_{i=9}^{16}  	&=& \frac{\rho^2}{2} \\
		    I_{i=17}^{20}  	&=&	\frac{\rho^3}{3} \\
		    I_{i=21}^{26}   &=&	\frac{\rho^4}{4} \\
		    I_{i=27}^{28}  	&=&	\frac{\rho^5}{5} \\
		    I_{i=29}^{31}  	&=& \frac{1-e^{-\gamma\rho^2}}{2\gamma}\\
		    I_{i=32}^{34}  	&=& \frac{1-\left(\rho^{2} \gamma + 1\right) e^{- \rho^{2} \gamma}}{2 \gamma^{2}}\\
    	\f}
     */
    double I(int i, double, double);

    /*!
     * H returns a multiplier in each term of the sum in P-3. This is used in
     * conjunction with C in the function Pp this represents the product
     * rho^n
     * - i=0 is the second additive in the formula in reynolds
     * - i=1 is the third ...
         \f{eqnarray*}
	 	  H_{i=0}^8 		&=& \rho^2 \\
	 	  H_{i=9}^{16}  	&=& \rho^3 \\
	 	  H_{i=17}^{20}  	&=&	\rho^4 \\
	 	  H_{i=21}^{26}   	&=&	\rho^5 \\
	 	  H_{i=27}^{28}  	&=&	\rho^6 \\
	 	  H_{i=29}^{31}  	&=& \rho^3 e^{-\gamma \rho^2}\\
	 	  H_{i=32}^{34}  	&=& \rho^5 e^{-\gamma \rho^2}\\

	 	\f}
     */
    double H(int i, double egrho);
};

}

#endif // ! TPX_HELIUM4_H
