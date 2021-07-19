//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolumeFraction.h"

registerMooseObject("gibbsApp", VolumeFraction);

template<>
InputParameters
validParams<VolumeFraction>()
{
  InputParameters params = validParams<ElementAverageValue>();
  return params;
}

VolumeFraction::VolumeFraction(const InputParameters & parameters)
 :ElementAverageValue(parameters)
{
}
 
Real
VolumeFraction::computeQpIntegral()
{                             
  return (_u[_qp]*_u[_qp]);
}
