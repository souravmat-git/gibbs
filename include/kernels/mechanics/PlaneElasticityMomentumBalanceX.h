//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "PlaneElasticityBase.h"

class PlaneElasticityMomentumBalanceX;

template <>
InputParameters validParams<PlaneElasticityMomentumBalanceX>();

/**
 * Kernel to implement the momentum balance
 * in x-direction for a single phase material
 * The variable that this kernel acts on is the displacement field in the x-direction.
 * Further, this variable is coupled to the displacement in the y-direction
 */
class PlaneElasticityMomentumBalanceX : public PlaneElasticityBase
{
public:
  PlaneElasticityMomentumBalanceX(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  //Displacement gradient in y-direction
  const VariableGradient & _grad_uy;
  unsigned int _uy_var;
  
  //Body force in X-direction
  const MaterialProperty<Real> & _fx;
 
};
