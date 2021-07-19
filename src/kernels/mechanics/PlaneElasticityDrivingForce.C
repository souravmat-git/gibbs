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

#include "PlaneElasticityDrivingForce.h"

registerMooseObject("gibbsApp", PlaneElasticityDrivingForce);

template<>
InputParameters
validParams<PlaneElasticityDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("disp_x", "Displacement in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  return params;
}

PlaneElasticityDrivingForce::PlaneElasticityDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  _grad_ux(coupledGradient("disp_x")),
  _ux_var(coupled("disp_x")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  _elastic_energy_alpha(getMaterialProperty<Real>("elastic_energy_alpha")),
  _elastic_energy_beta(getMaterialProperty<Real>("elastic_energy_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _L(getMaterialProperty<Real>("mob_name"))
{
}

Real
PlaneElasticityDrivingForce::computeQpResidual(){
  //Difference in elastic strain energy density is the driving force of PT
  //Residual: N_{i} * h_{\phi} (f_el_beta - f_el_alpha) = 0 ;
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * (_elastic_energy_beta[_qp] - _elastic_energy_alpha[_qp]);
}

Real
PlaneElasticityDrivingForce::computeQpJacobian(){
  // This will be non-zero since the 
  // free energies of each phase is dependent strain
  return  _test[_i][_qp] * _L[_qp] * _d2h[_qp] * (_elastic_energy_beta[_qp] - _elastic_energy_alpha[_qp]) *_phi[_j][_qp];
}

Real
PlaneElasticityDrivingForce::computeQpOffDiagJacobian(unsigned int jvar){
  if (jvar == _ux_var)
  {
    return _test[_i][_qp] *_L[_qp] * _dh[_qp] *(_sx_beta[_qp] - _sx_alpha[_qp])* _grad_phi[_j][_qp](0);  
  }
  else
    return 0.0;
}
