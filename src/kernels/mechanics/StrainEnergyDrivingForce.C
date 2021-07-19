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

#include "StrainEnergyDrivingForce.h"

registerMooseObject("gibbsApp", StrainEnergyDrivingForce);

template<>
InputParameters
validParams<StrainEnergyDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");
  params.addRequiredCoupledVar("disp_y", "Displacement in the y-direction");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addParam<MaterialPropertyName>("nd_factor",1.0, "RT/(Vm*barrier_height)");
  return params;
}

StrainEnergyDrivingForce::StrainEnergyDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  //_grad_ux(coupledGradient("disp_x")),
  _ux_var(coupled("disp_x")),
  //_grad_uy(coupledGradient("disp_y")),
  _uy_var(coupled("disp_y")),
  //Interpolation function
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  //Requisite materila properties
  _fel_alpha(getMaterialProperty<Real>("fel_alpha")),
  _fel_beta(getMaterialProperty<Real>("fel_beta")),
  //Normal stress in the x-direction
  _sxx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sxx_beta(getMaterialProperty<Real>("sx_beta")),
  //Normal stress in the y-direction
  _syy_alpha(getMaterialProperty<Real>("sy_alpha")),
  _syy_beta(getMaterialProperty<Real>("sy_beta")),
  //Shear stress in the xy-direction
  _sxy_alpha(getMaterialProperty<Real>("sxy_alpha")),
  _sxy_beta(getMaterialProperty<Real>("sxy_beta")),
  //Mobility and non-dimensional factor
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
StrainEnergyDrivingForce::computeQpResidual()
{
  //Difference in elastic strain energy density is the driving force of PT
  //Residual: N_{i} * h_{\phi} (f_el_beta - f_el_alpha) = 0 ;
  return  _nd_factor[_qp] * _test[_i][_qp] *_L[_qp] * _dh[_qp] * (_fel_beta[_qp] - _fel_alpha[_qp]);

}

Real
StrainEnergyDrivingForce::computeQpJacobian()
{
  // This will be non-zero since the 
  // free energies of each phase is dependent strain
  return  _nd_factor[_qp] * _test[_i][_qp] *_L[_qp] * _d2h[_qp] * (_fel_beta[_qp] - _fel_alpha[_qp]) *_phi[_j][_qp];
}

Real
StrainEnergyDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ux_var)
  {
    return _nd_factor[_qp] * _test[_i][_qp] * _L[_qp] * _dh[_qp] * (_sxx_beta[_qp] - _sxx_alpha[_qp]) * _grad_phi[_j][_qp](0)
          +_nd_factor[_qp] * _test[_i][_qp] * _L[_qp] * _dh[_qp] * (_sxy_beta[_qp] - _sxy_alpha[_qp]) * _grad_phi[_j][_qp](1);  
  }
  else if (jvar == _uy_var)
  {
    return _nd_factor[_qp] * _test[_i][_qp] * _L[_qp] * _dh[_qp] * (_sxy_beta[_qp] - _sxy_alpha[_qp]) * _grad_phi[_j][_qp](0)
          +_nd_factor[_qp] * _test[_i][_qp] * _L[_qp] * _dh[_qp] * (_syy_beta[_qp] - _syy_alpha[_qp]) * _grad_phi[_j][_qp](1);  
  }
  else
    return 0.0;
}
