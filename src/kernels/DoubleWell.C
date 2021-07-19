//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This was written by S.Chatterjee
//f(phi) = phi^(2) * (1-phi)^(2)

#include "DoubleWell.h"

registerMooseObject("gibbsApp", DoubleWell);

template <>
InputParameters
validParams<DoubleWell>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the barrier term"
                              "Eqn: H*dg/dphi");
  params.addRequiredParam<MaterialPropertyName>("barrier_height", "Height of the double well potential");
  params.addRequiredParam<MaterialPropertyName>("mob_name", "The mobility used with the kernel");

  return params;
}

DoubleWell::DoubleWell(const InputParameters & parameters)
  : Kernel(parameters),
    _BH(getMaterialProperty<Real>("barrier_height")),
    _L(getMaterialProperty<Real>("mob_name"))
{
}

Real
DoubleWell::computeQpResidual()
{
  return _L[_qp] * _test[_i][_qp] * _BH[_qp] * 2 *_u[_qp] * (1-_u[_qp]) * (1- 2 * _u[_qp]);
}

Real
DoubleWell::computeQpJacobian()
{
  return  _L[_qp] * _test[_i][_qp] *_BH[_qp] *_phi[_j][_qp]* 2 * (6 * _u[_qp] * _u[_qp] - 6 * _u[_qp] + 1);
}

Real
DoubleWell::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0.0;
}
