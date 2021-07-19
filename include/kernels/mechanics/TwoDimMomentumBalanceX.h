//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Kernel.h"

class TwoDimMomentumBalanceX;

template <>
InputParameters validParams<TwoDimMomentumBalanceX>();

/**
 * Kernel to implement the momentum balance
 * in x-direction for a single phase material
 * Note that this code is not dimension agnostic
 * and depends on the load
 */
class TwoDimMomentumBalanceX : public Kernel
{
public:
  TwoDimMomentumBalanceX(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
    
  //Youngs Modulus of a given phase
  const MaterialProperty<Real> & _E;
  
  //Body force in X-direction
  const MaterialProperty<Real> & _bx;
  
  //For non-dimensionalization we will need two parameters
  const MaterialProperty<Real> & _lc;
   
};
