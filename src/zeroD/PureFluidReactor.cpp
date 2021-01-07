//! @file PureFluidReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/PureFluidReactor.h"

using namespace std;

namespace Cantera
{

void PureFluidReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "PureFluid") {
        throw CanteraError("PureFluidReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void PureFluidReactor::updateState(doublereal* y)
{
	// The components of y are [0] the total mass, [1] the total volume,
	// [2] the total internal energy, [3...K+3] are the mass fractions of each
	// species, and [K+3...] are the coverages of surface species on each wall.
	m_mass = y[0];
	m_vol = y[1];
	m_thermo->setMassFractions_NoNorm(y+3);

	if (m_energy) {
		m_thermo->setState_UV(y[2]/m_mass, m_vol/m_mass);
	} else {
		m_thermo->setDensity(m_mass/m_vol);
	}

	updateSurfaceState(y + m_nsp + 3);

	// save parameters needed by other connected reactors
	m_enthalpy = m_thermo->enthalpy_mass();
	m_pressure = m_thermo->pressure();
	m_intEnergy = m_thermo->intEnergy_mass();
	m_thermo->saveState(m_state);
}

}
