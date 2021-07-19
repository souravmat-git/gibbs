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

#include "KHSDrivingForce.h"
registerMooseObject("gibbsApp", KHSDrivingForce);

template<>
InputParameters
validParams<KHSDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor","Non-dimensional factor");
  return params;
}

KHSDrivingForce::KHSDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _ux_var(coupled("disp_x")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  _elastic_energy_alpha(getMaterialProperty<Real>("elastic_energy_alpha")),
  _elastic_energy_beta(getMaterialProperty<Real>("elastic_energy_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
KHSDrivingForce::computeQpResidual()
{
  //Difference in elastic strain energy density is the driving force of PT
  //Residual: N_{i} * h_{\phi} (f_el_beta - f_el_alpha) = 0 ;
  
  Real diff_el = (_elastic_energy_beta[_qp] - _elastic_energy_alpha[_qp]); 
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * diff_el * _nd_factor[_qp];

}

Real
KHSDrivingForce::computeQpJacobian()
{
  // This will be non-zero since the 
  // free energies of each phase is dependent strain
  Real diff_el = (_elastic_energy_beta[_qp] - _elastic_energy_alpha[_qp]); 
  return  _test[_i][_qp] * _L[_qp] * _d2h[_qp] * diff_el *_phi[_j][_qp] * _nd_factor[_qp];
}

Real
KHSDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ux_var){ 
    Real diff_sx = (_sx_beta[_qp] - _sx_alpha[_qp]); 
    return _test[_i][_qp] *_L[_qp] * _dh[_qp] * diff_sx * _grad_phi[_j][_qp](0)*_nd_factor[_qp];  
  }
  else
    return 0.0;
}
