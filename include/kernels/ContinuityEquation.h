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

class ContinuityEquation;

template <>
InputParameters validParams<ContinuityEquation>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 */
class ContinuityEquation : public Kernel
{
public:
  ContinuityEquation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


private:

  //Declare constant member function
  Real thermodynamic_factor() const;
  Real L_BB_interp() const;
  Real dL_BB_xB_interp() const;
  
  const VariableValue & _eta;
  unsigned int _eta_var;
  
  const VariableGradient &  _xB;
  unsigned int _xB_var;
  
  //Thermodynamic factor
  const MaterialProperty<Real> & _B_tf_alpha;
  const MaterialProperty<Real> & _B_tf_beta;
  
  //Interpolation functions
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  /// Isotropic diffusion mobility dependent of concentration
  const MaterialProperty<Real> & _L_BB_alpha;
  const MaterialProperty<Real> & _L_BB_beta;
  
  const MaterialProperty<Real> & _dL_BB_xB_alpha;
  const MaterialProperty<Real> & _dL_BB_xB_beta;
};
//#endif // CONTINUITYEQUATION_H
