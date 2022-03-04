//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*Written by S.Chatterjee*/

#include "PartialRankOneHomogenization.h"
registerMooseObject("gibbsApp", PartialRankOneHomogenization);

//template <>
InputParameters
PartialRankOneHomogenization::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
                                            "Name of strain jump material");
  params.addRequiredCoupledVar("eta", "phase-field variable");
  params.addClassDescription("Computes the strain jump [e] tensor");
  return params;
}

PartialRankOneHomogenization::PartialRankOneHomogenization(const InputParameters & parameters)
  : Material(parameters),
   //Magnitude of Jump and its derivative wrt h
   _a(getMaterialProperty<RealVectorValue>("a")),
   _da_dphi(getMaterialProperty<RealVectorValue>("da_dphi")),
    //Gradient of phase-field variable
   _grad_eta(coupledGradient("eta")),
   //Unit normal
   //_n(getMaterialProperty<RealVectorValue>("n")),
   //The (compatible) strain jump tensor will be calculated by this class
   _strain_jump(declareProperty<RankTwoTensor>
                          (getParam<MaterialPropertyName>("strain_jump_name"))),
   _dstrainjump_dphi(declareProperty<RankTwoTensor>("dstrainjump_dphi"))
{
}

void
PartialRankOneHomogenization::computeQpProperties(){

  const Real _norm_sq = _grad_eta[_qp].norm_sq();

  if (_norm_sq < std::pow(libMesh::TOLERANCE, 3.0)){

     _strain_jump[_qp].zero();
     _dstrainjump_dphi[_qp].zero();

  }
  else {
   //For partial rank-one homogenization
   //0. Declare the unit normal vector
   //1. Calculate the jump in displacement gradient
   //2. And then the symmetric part of that tensor

   RealVectorValue _n = - _grad_eta[_qp]/std::sqrt(_norm_sq);

   RankTwoTensor _jump_in_disp_grad_tensor, _deriv_jump_in_disp_grad_tensor;

  _jump_in_disp_grad_tensor.vectorOuterProduct(_a[_qp], _n);
  _deriv_jump_in_disp_grad_tensor.vectorOuterProduct(_da_dphi[_qp], _n);

  _strain_jump[_qp] = (_jump_in_disp_grad_tensor
                    +  _jump_in_disp_grad_tensor.transpose()) / 2.0;

  _dstrainjump_dphi[_qp] = (_deriv_jump_in_disp_grad_tensor
                        +  _deriv_jump_in_disp_grad_tensor.transpose()) / 2.0;
  }
}
