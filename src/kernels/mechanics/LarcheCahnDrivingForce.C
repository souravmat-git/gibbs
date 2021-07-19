//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html 

//* This kernel implements the driving force 
//* which is the difference in elastic strain energy between the two phases

#include "LarcheCahnDrivingForce.h"

registerMooseObject("gibbsApp", LarcheCahnDrivingForce);

template<>
InputParameters
validParams<LarcheCahnDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("stress_x", "Stress in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-Dimensional factor");
  return params;
}

LarcheCahnDrivingForce::LarcheCahnDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  _stress_x(coupledValue("stress_x")),
  _stress_x_var(coupled("stress_x")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  _comp_energy_alpha(getMaterialProperty<Real>("comp_energy_alpha")),
  _comp_energy_beta(getMaterialProperty<Real>("comp_energy_beta")),
  _ex_alpha(getMaterialProperty<Real>("ex_alpha")),
  _ex_beta(getMaterialProperty<Real>("ex_beta")),
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
LarcheCahnDrivingForce::computeQpResidual()
{
  //Difference in elastic strain energy density is the driving force of PT
  //Residual: N_{i} * h_{\phi} (f_el_beta - f_el_alpha) = 0 ;
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * _nd_factor[_qp]*(_comp_energy_beta[_qp] - _comp_energy_alpha[_qp]);
}

Real
LarcheCahnDrivingForce::computeQpJacobian()
{
  // This will be non-zero since the 
  // free energies of each phase is dependent strain
  return  _test[_i][_qp] * _L[_qp] * _d2h[_qp] * _nd_factor[_qp]*(_comp_energy_beta[_qp] - _comp_energy_alpha[_qp]) *_phi[_j][_qp];
}

Real
LarcheCahnDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _stress_x_var)
  {
    return _test[_i][_qp] *_L[_qp] * _dh[_qp] * _nd_factor[_qp]*(-_ex_beta[_qp] + _ex_alpha[_qp]) * _phi[_j][_qp];  
  }
  else
    return 0.0;
}
