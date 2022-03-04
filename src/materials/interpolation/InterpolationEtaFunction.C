//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This Material class supplies the interpolation function,
//* and its first and second derivatives

#include "InterpolationEtaFunction.h"
#include "Function.h"
registerMooseObject("gibbsApp", InterpolationEtaFunction);

//template <>
InputParameters
InterpolationEtaFunction::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<FunctionName>("eta", "The coupled variable to the material property"
                       "function");
  return params;
}

InterpolationEtaFunction::InterpolationEtaFunction(const InputParameters & parameters)
  : Material(parameters),
    _h(declareProperty<Real>("h")),
    _dh(declareProperty<Real>("dh")),
    _d2h(declareProperty<Real>("d2h")),
    //_g(declareProperty<Real>("g")),
    _eta(getFunction("eta"))
{
}

void
InterpolationEtaFunction::computeQpProperties()
{
   Real _phi = _eta.value(_t, _q_point[_qp]);

  _h[_qp] = std::pow(_phi,3.0) * (6.0 * std::pow(_phi,2.0) - 15.0 * _phi + 10);

  _dh[_qp] = 30.0 * std::pow(_phi,2.0) * std::pow((1.0-_phi),2.0);

  _d2h[_qp] = 60.0 * _phi * (1.0-_phi) * (1.0 - 2.0 * _phi);

}
