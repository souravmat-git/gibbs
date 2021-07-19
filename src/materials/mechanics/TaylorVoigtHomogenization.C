//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "TaylorVoigtHomogenization.h"
registerMooseObject("gibbsApp",TaylorVoigtHomogenization);

template <>
InputParameters
validParams<TaylorVoigtHomogenization>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  params.addClassDescription("Computes the strain jump [e] tensor");
  return params;
}

TaylorVoigtHomogenization::TaylorVoigtHomogenization(const InputParameters & parameters)
  : Material(parameters),
   //The (compatible) strain jump tensor will be calculated by this class
   _strain_jump(declareProperty<RankTwoTensor>
                          (getParam<MaterialPropertyName>("strain_jump_name"))),
   _dstrainjump_dphi(declareProperty<RankTwoTensor>("dstrainjump_dphi"))
{
} 

void
TaylorVoigtHomogenization::computeQpProperties(){  
   
  _strain_jump[_qp].zero();
  _dstrainjump_dphi[_qp].zero();

}

