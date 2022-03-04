//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputePhaseStress.h"
registerMooseObject("gibbsApp", ComputePhaseStress);

//template <>
InputParameters
ComputePhaseStress::validParams()
{
  InputParameters params = ComputeStressBase::validParams();
  return params;
}

ComputePhaseStress::ComputePhaseStress(const InputParameters & parameters)
  : ComputeStressBase(parameters),
   //elasticity tensor to calculate stress
   _stiffness(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor"))
{
}

void
ComputePhaseStress::computeQpStress(){

   _stress[_qp] = _stiffness[_qp] * _mechanical_strain[_qp];
   _elastic_strain[_qp] = _mechanical_strain[_qp];

}
