//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceEnergy.h"

registerMooseObject("gibbsApp", InterfaceEnergy);

template<>
InputParameters
validParams<InterfaceEnergy>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("phase_beta", "Phase field variable beta, assumed to be eta=1");
  params.addCoupledVar("phase_alpha", 0.0, "Phase field variable alpha, assumed to be eta=0");
  params.addParam<MaterialPropertyName>("barrier_height", 1.0, "Height of double well");
  return params;
}

InterfaceEnergy::InterfaceEnergy(const InputParameters & parameters)
  :AuxKernel(parameters),
  _phase_alpha(coupledValue("phase_alpha")),
  _grad_phase_alpha(coupledGradient("phase_alpha")),
  _phase_beta(coupledValue("phase_beta")),
  _grad_phase_beta(coupledGradient("phase_beta")),
  _kappa(getMaterialProperty<Real>("kappa")),
  _BH(getMaterialProperty<Real>("barrier_height"))
{
}
 
Real
InterfaceEnergy::computeValue()
{
  const Real gradient_energy = 0.5* _kappa[_qp] * _grad_phase_alpha[_qp](0) * _grad_phase_alpha[_qp](0)
                            +  0.5* _kappa[_qp] * _grad_phase_beta[_qp](0) * _grad_phase_beta[_qp](0) ;
  
  const Real doublewell = _BH[_qp]* _phase_alpha[_qp]* _phase_alpha[_qp]* std::pow((1.0 - _phase_alpha[_qp]),2.0)+
                          _BH[_qp]* _phase_beta[_qp]* _phase_beta[_qp]* std::pow((1.0 - _phase_beta[_qp]),2.0);
                             
  return (gradient_energy + doublewell);
}
