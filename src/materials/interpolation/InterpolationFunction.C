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

#include "InterpolationFunction.h"
registerMooseObject("gibbsApp", InterpolationFunction);

//template <>
InputParameters
InterpolationFunction::validParams()
{
  InputParameters params = Material::validParams();
  params.addCoupledVar("eta", "The coupled variable to the material property"
                       "function");
  return params;
}

InterpolationFunction::InterpolationFunction(const InputParameters & parameters)
  : Material(parameters),
    // Declare that this material is going to provide a Real
    // valued properties named "_h,_dh,_d2h" that Kernels can use.
    _h(declareProperty<Real>("h")),
    _dh(declareProperty<Real>("dh")),
    _d2h(declareProperty<Real>("d2h")),
    //_g(declareProperty<Real>("g")),
    _eta(coupledValue("eta"))
{
}


void
InterpolationFunction::computeQpProperties()
{
  //The interpolation function is of the form
  _h[_qp] = std::pow(_eta[_qp],3.0) * (6.0 * std::pow(_eta[_qp],2.0) - 15.0 * _eta[_qp] + 10);

  // First derivative of the interpolation function
  _dh[_qp] = 30.0 * std::pow(_eta[_qp],2.0) * std::pow((1.0-_eta[_qp]),2.0);

  // Second derivative of the interpolation function
  _d2h[_qp] = 60.0 * _eta[_qp] * (1.0-_eta[_qp]) * (1.0 - 2.0 * _eta[_qp]);

  //_g[_qp] = _eta[_qp] * _eta[_qp] * std::pow((1.0-_eta[_qp]),2.0);
}
