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

#include "QDrivingForce.h"
registerMooseObject("gibbsApp", QDrivingForce);

template<>
InputParameters
validParams<QDrivingForce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: dh*[f_beta_el -f_alpha_el] = 0");
  params.addRequiredCoupledVar("ux", "Displacecment in the x-direction");
  params.addRequiredParam<MaterialPropertyName>("a", "Jump in strain");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "phase field mobility");
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

QDrivingForce::QDrivingForce(const InputParameters & parameters)
  :Kernel(parameters),
  _ux_var(coupled("ux")),
  //Interpolation function
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  _d2h(getMaterialProperty<Real>("d2h")),
  //Requisite material properties
  _fel_alpha(getMaterialProperty<Real>("fel_alpha")),
  _fel_beta(getMaterialProperty<Real>("fel_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  //Strain jump is required for this kernel
  _a(getMaterialProperty<Real>(getParam<MaterialPropertyName>("a"))),
  //modulus values
  _mat_const_alpha(getMaterialProperty<Real>("mat_const_alpha")),
  _mat_const_beta(getMaterialProperty<Real>("mat_const_beta")),
  _L(getMaterialProperty<Real>("mob_name")),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

Real
QDrivingForce::fdf() const{
  return (_fel_beta[_qp] - _fel_alpha[_qp]) - 
             _a[_qp]*(_sx_beta[_qp]* _h[_qp] + _sx_alpha[_qp]* (1.0- _h[_qp]));
}

Real
QDrivingForce::computeQpResidual(){ 
  return _test[_i][_qp] *_L[_qp] * _dh[_qp] * fdf() * _nd_factor[_qp];

}

Real
QDrivingForce::computeQpJacobian(){
  return _test[_i][_qp]* _L[_qp] * _d2h[_qp] * fdf() * _nd_factor[_qp]* _phi[_j][_qp];
}

Real
QDrivingForce::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _ux_var){
    
     Real _dfdf_de = _a[_qp] * _sx_beta[_qp]*(_mat_const_beta[_qp])/
                           ((1.0 -_h[_qp])* _mat_const_beta[_qp] + _h[_qp] * _mat_const_alpha[_qp]);
     
     return _test[_i][_qp] *_L[_qp] * _dh[_qp] * _dfdf_de * _nd_factor[_qp] * _grad_phi[_j][_qp](0);
  }
  else
    return 0;   
}
