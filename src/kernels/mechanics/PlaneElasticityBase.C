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

#include "PlaneElasticityBase.h"

registerMooseObject("gibbsApp", PlaneElasticityBase);

template <>
InputParameters
validParams<PlaneElasticityBase>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("C11", "Stiffness constants");
  params.addRequiredParam<MaterialPropertyName>("C12", "Stiffness constants");
  params.addRequiredParam<MaterialPropertyName>("C22", "Stiffness constants");
  params.addRequiredParam<MaterialPropertyName>("C66", "Stiffness constants");
  params.addClassDescription("This kernel is the base kernel for plane elasticity problem");

  return params;
}

PlaneElasticityBase::PlaneElasticityBase(const InputParameters & parameters)
  :Kernel(parameters),
  _C11(getMaterialProperty<Real>(getParam<MaterialPropertyName>("C11"))),
  _C12(getMaterialProperty<Real>(getParam<MaterialPropertyName>("C12"))),
  _C22(getMaterialProperty<Real>(getParam<MaterialPropertyName>("C22"))),
  _C66(getMaterialProperty<Real>(getParam<MaterialPropertyName>("C66")))
{
}

Real
PlaneElasticityBase::computeQpResidual(){
  return 0;
}

Real
PlaneElasticityBase::computeQpJacobian(){
  //Derivative with respect to the non-linear variable
  return 0;
}

Real
PlaneElasticityBase::computeQpOffDiagJacobian(unsigned int /*jvar*/){ 
    return 0;
}
