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

    virtual void updateState(doublereal* y);

};

}

#endif
