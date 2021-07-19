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

#include "LagrangeMultiplierPF.h"

registerMooseObject("gibbsApp", LagrangeMultiplierPF);

template <>
InputParameters
validParams<LagrangeMultiplierPF>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Lagrange Multiplier for phase field");
  params.addRequiredCoupledVar("lambda", "lagrange multiplier");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  return params;
}

LagrangeMultiplierPF::LagrangeMultiplierPF(const InputParameters & parameters)
  : Kernel(parameters),
   _lambda(coupledValue("lambda")),
   _lambda_var(coupled("lambda")),
   _L(getMaterialProperty<Real>("mob_name"))
{
}

Real
LagrangeMultiplierPF::computeQpResidual()
{
  return -(_L[_qp] * _test[_i][_qp] * _lambda[_qp]);
}

Real
LagrangeMultiplierPF::computeQpJacobian()
{
  return  0.0;
}

Real
LagrangeMultiplierPF::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _lambda_var)
    return -(_L[_qp] * _test[_i][_qp] * _phi[_j][_qp]);
  else
    return 0;
  
}
