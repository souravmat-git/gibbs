//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "LamesParameters.h"
registerMooseObject("gibbsApp", LamesParameters);

//template <>
InputParameters
LamesParameters::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::string>("base_name","The phase name");
  return params;
}

LamesParameters::LamesParameters(const InputParameters & parameters)
  : Material(parameters),
    _phase_name(getParam<std::string>("base_name")),
   //Stiffness of alpha and beta phase
   _stiffness(getMaterialProperty<RankFourTensor>(_phase_name + "_elasticity_tensor")),
   //Compute the following properties
   _lambda(declareProperty<Real>(_phase_name + "_lambda")),
   _mu(declareProperty<Real>(_phase_name + "_mu"))
{
}

void
LamesParameters::computeQpProperties()
{
    //Initialize the Lames constant //C_1122 = C12
   _lambda[_qp] = _stiffness[_qp](0,0,1,1);

   //mu = C_2323 = C_44 (Voigt Notation)
   _mu[_qp]     = _stiffness[_qp](1,2,1,2);
}
