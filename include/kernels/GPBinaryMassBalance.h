//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GPBINARYMASSBALANCE_H
#define GPBINARYMASSBALANCE_H

#include "Kernel.h"

class GPBinaryMassBalance;

template <>
InputParameters validParams<GPBinaryMassBalance>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * variable on which this kernel operates: mu
 */
class GPBinaryMassBalance: public Kernel
{
public:
  GPBinaryMassBalance(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:  

  //Declare constant memeber functon
  Real L_BB_interp() const;
  Real dL_BB_muB_interp() const;
 
  const VariableValue & _eta;
  unsigned int _eta_var;

  //_h is a reference variable that indicates the interpolation function
  const MaterialProperty<Real>& _h;
  // _dh is a reference variable that indicates the first derivative of the 
  // interpolation function
  const MaterialProperty<Real>& _dh; 
  
  /// Isotropic diffusion mobility dependent on diffusion potential
  const MaterialProperty<Real> & _L_BB_beta;
  
  const MaterialProperty<Real> & _L_BB_alpha;
  
  const MaterialProperty<Real> & _dL_BB_muB_beta;
  
  const MaterialProperty<Real> & _dL_BB_muB_alpha;
};

#endif // GPBINARYMASSBALANCE_H
