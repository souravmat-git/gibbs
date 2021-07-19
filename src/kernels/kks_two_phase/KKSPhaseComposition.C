//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This code modifies the KKSPhaseConcentration in MOOSE
//* by including the hand coded interpolation function
//* instead of using the Material property- S.Chatterjee

#include "KKSPhaseComposition.h"

registerMooseObject("gibbsApp", KKSPhaseComposition);

template <>
InputParameters
validParams<KKSPhaseComposition>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (1-h(eta))*c_alpha + h(eta)*c_beta - c = 0."
                             "non-linear variable of this kernel is c_beta.");
  params.addRequiredCoupledVar("xB_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("mole_fraction", "Real concentration");
  params.addRequiredCoupledVar("eta", "To distinguish phases");
  return params;
}

KKSPhaseComposition::KKSPhaseComposition(const InputParameters & parameters)
  :  Kernel(parameters),
    _xB_alpha(coupledValue("xB_alpha")),
    _xB_alpha_var(coupled("xB_alpha")),
    _mole_fraction(coupledValue("mole_fraction")),
    _mole_fraction_var(coupled("mole_fraction")),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _h(getMaterialProperty<Real>("h")),
    _dh(getMaterialProperty<Real>("dh"))
{
}

Real
KKSPhaseComposition::computeQpResidual()
{
  // R: w*((1-h(eta))*ca + h(eta)*cb - c)
  
  return _test[_i][_qp] * (_h[_qp] *_u[_qp] + (1-_h[_qp]) * _xB_alpha[_qp] 
         - _mole_fraction[_qp]);
}

Real
KKSPhaseComposition::computeQpJacobian()
{
 
  return _test[_i][_qp] * _h[_qp] * _phi[_j][_qp];
}

Real
KKSPhaseComposition::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _xB_alpha_var)
  {
    return _test[_i][_qp] * (1.0 - _h[_qp]) * _phi[_j][_qp];
  }
  else if (jvar == _mole_fraction_var)
  {
    return -_test[_i][_qp] * _phi[_j][_qp];
  }
  else if (jvar == _eta_var)
  {
    return _test[_i][_qp] * (_u[_qp] - _xB_alpha[_qp]) * _dh[_qp] * _phi[_j][_qp];
  }  
  else
    return 0.0;
}
