//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "QElasticDrivingForce.h"
registerMooseObject("gibbsApp", QElasticDrivingForce);

template<>
InputParameters
validParams<QElasticDrivingForce>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<MaterialPropertyName>("a", "Jump in strain");
  return params;
}

QElasticDrivingForce::QElasticDrivingForce(const InputParameters & parameters)
  :AuxKernel(parameters),
  //Interpolation function
  _h(getMaterialProperty<Real>("h")),
  _dh(getMaterialProperty<Real>("dh")),
  //Requisite material properties
  _fel_alpha(getMaterialProperty<Real>("fel_alpha")),
  _fel_beta(getMaterialProperty<Real>("fel_beta")),
  _sx_alpha(getMaterialProperty<Real>("sx_alpha")),
  _sx_beta(getMaterialProperty<Real>("sx_beta")),
  _a(getMaterialProperty<Real>(getParam<MaterialPropertyName>("a"))),
  _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}
 
Real
QElasticDrivingForce::computeValue()
{
  
  return _dh[_qp] * _nd_factor[_qp] * ( (_fel_beta[_qp] - _fel_alpha[_qp]) - 
                                        _a[_qp]*(_sx_beta[_qp]* _h[_qp] + _sx_alpha[_qp]* (1.0- _h[_qp])));
}
