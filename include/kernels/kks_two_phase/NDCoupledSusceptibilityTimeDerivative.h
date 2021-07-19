//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CoupledTimeDerivative.h"
#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"
#include "CoupledSusceptibilityTimeDerivative.h"

// Forward Declaration
class NDCoupledSusceptibilityTimeDerivative;

template <>
InputParameters validParams<NDCoupledSusceptibilityTimeDerivative>();

/**
 * This calculates a modified coupled time derivative that multiplies the time derivative of a
 * coupled variable by a function of the variables
 */
class NDCoupledSusceptibilityTimeDerivative
    : public CoupledSusceptibilityTimeDerivative
{
public:
  NDCoupledSusceptibilityTimeDerivative(const InputParameters & parameters);


protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
private:
  
 //phase comp_alpha is a function of diffusion potential
  const MaterialProperty<Real> & _xB_alpha;
  
  //Phase comp_beta is a function of diffusion potential
  const MaterialProperty<Real> & _xB_beta;
  
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> &  _inv_B_tf_beta;
  
  // _dh is a reference variable that indicates the first derivative of the 
  // interpolation function
  const MaterialProperty<Real>& _dh;
  const MaterialProperty<Real>& _d2h;  
};

