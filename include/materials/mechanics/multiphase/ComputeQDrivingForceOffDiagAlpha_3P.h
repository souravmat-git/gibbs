//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class ComputeQDrivingForceOffDiagAlpha_3P;

//MOOSe includes
#include "ComputeQDrivingForceOffDiagBase_3P.h"

template <>
InputParameters validParams<ComputeQDrivingForceOffDiagAlpha_3P>();

class ComputeQDrivingForceOffDiagAlpha_3P : public ComputeQDrivingForceOffDiagBase_3P
{
public:
  ComputeQDrivingForceOffDiagAlpha_3P(const InputParameters & parameters);
  //This class calculates the off diag term corresponding to the alpha phase

protected:
  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Derivative interpolation
    const MaterialProperty<Real> & _dhalpha_dphialpha;
    const MaterialProperty<Real> & _dhbeta_dphialpha;
    const MaterialProperty<Real> & _dhgamma_dphialpha;

    MaterialProperty<RankTwoTensor> & _d2Fdcdstrain_alpha;
};
