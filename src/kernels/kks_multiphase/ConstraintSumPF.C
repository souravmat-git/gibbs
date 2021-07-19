//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*In this code we explicitly code the free energy
//*The form of the free energy is f(eta) phi^(4)/4 - phi^(2)/2
//The non-linear variable that this kernel operates on is eta

#include "ConstraintSumPF.h"

registerMooseObject("gibbsApp", ConstraintSumPF);

template <>
InputParameters
validParams<ConstraintSumPF>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Sum of phase fields = 1");
  params.addRequiredCoupledVar("phase_alpha", "Variable representing alpha phase");
  params.addRequiredCoupledVar("phase_beta", "Variable representing beta phase");
  params.addCoupledVar("phase_gamma", 0.0, "Variable representing beta phase");
  return params;
}

ConstraintSumPF::ConstraintSumPF(const InputParameters & parameters)
  : Kernel(parameters),
 _phase_alpha(coupledValue("phase_alpha")),
 _phase_alpha_var(coupled("phase_alpha")),
 _phase_beta(coupledValue("phase_beta")),
 _phase_beta_var(coupled("phase_beta")),
 _phase_gamma(coupledValue("phase_gamma")),
 _phase_gamma_var(coupled("phase_gamma"))
{
}

Real
ConstraintSumPF::computeQpResidual()
{
  return (_test[_i][_qp]  * (1.0 - (_phase_alpha[_qp] + _phase_beta[_qp] + _phase_gamma[_qp])) );
}

Real
ConstraintSumPF::computeQpJacobian()
{ 
  return  (0.0);
}

Real
ConstraintSumPF::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _phase_alpha_var)
  { 
    return -(_test[_i][_qp]  * _phi[_j][_qp]);
  }
  else if (jvar == _phase_beta_var)
  {
    return -( _test[_i][_qp] * _phi[_j][_qp]);
  }
  else if (jvar == _phase_gamma_var)
  {
    return -( _test[_i][_qp] * _phi[_j][_qp]);
  }  
  else 
    return 0; 
}
