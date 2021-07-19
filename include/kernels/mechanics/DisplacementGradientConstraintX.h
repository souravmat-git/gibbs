//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//Include dependencies
#include "ObtainUnitNormalBase.h"

// Forward Declarations
class DisplacementGradientConstraintX;

template <>
InputParameters validParams<DisplacementGradientConstraintX>();


class DisplacementGradientConstraintX : public ObtainUnitNormalBase
{
public:
  DisplacementGradientConstraintX(const InputParameters & parameters);

protected:  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  Real axx() const;
  
  const VariableValue & _eta;
  unsigned int _eta_var;

  const VariableValue  & _ex_alpha;
  unsigned int _ex_alpha_var;

  const VariableGradient & _grad_ux;
  unsigned int _ux_var;
 
  const MaterialProperty<Real> & _eT_alpha;
  const MaterialProperty<Real> & _eT_beta;
  
   const MaterialProperty<Real> & _h;
   const MaterialProperty<Real> & _dh;
};
