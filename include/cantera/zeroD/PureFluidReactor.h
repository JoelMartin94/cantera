//! @file PureFluidReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PUREFLUIDREACTOR_H
#define CT_PUREFLUIDREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class PureFluidReactor is a class for stirred reactors that is specifically
 * optimized for pure fluids.
 */
class PureFluidReactor : public Reactor
{
public:
    PureFluidReactor() {}

    virtual int type() const {
        return PureFluidReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);

protected:
    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
