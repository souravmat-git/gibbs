//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableUnitNormal.h"
registerMooseObject("gibbsApp", VariableUnitNormal);

//template<>
InputParameters
VariableUnitNormal::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Calculates the projection tensor from the PFV");
  params.addParam<MaterialPropertyName>("unit_normal_name","n", "Unit normal name");
  //params.addParam<Real>("tol", 1e-20, "Tolerance for the gradient");
  params.addRequiredCoupledVar("eta", "Phase-field variable");
  params.addRequiredParam<RealVectorValue>("n", "A value in case eta gradient is zero"
                                         "in the bulk phases");
  return params;
}

VariableUnitNormal::VariableUnitNormal(const InputParameters & parameters)
 : Material(parameters),
  _grad_eta(coupledGradient("eta")),
  _n(getParam<RealVectorValue>("n")),
  _unit_normal_name(getParam<MaterialPropertyName>("unit_normal_name")),
  //_tol(getParam<Real>("tol")),
  _n_val(declareProperty<RealVectorValue>(_unit_normal_name)),
  _dn_dgradphi(declareProperty<RankTwoTensor>("dn_dgradphi"))
{
}

void
VariableUnitNormal::computeQpProperties()
{
  const Real _norm_sq   =  _grad_eta[_qp].norm_sq();
  const Real _norm_cube =  std::pow(_grad_eta[_qp].norm(),3.0);

  const RankTwoTensor I(RankTwoTensor::initIdentity);
  RankTwoTensor OG;

  OG.vectorOuterProduct(_grad_eta[_qp], _grad_eta[_qp]);

  if (_norm_sq < std::pow(libMesh::TOLERANCE,3.0)){
   //Since unit normal is constant its derivative is zero
   _n_val[_qp] = _n;
   _dn_dgradphi[_qp] = 0;

   //std::cout << _n_val[_qp].norm_sq() << std::endl;
   //std::cout << libMesh::TOLERANCE <<  std::endl;

  }
  else{
   //The unit normal is defined such that it points from the
   //beta phase (where phi = 1) to alpha phase (where phi = 0)
   _n_val[_qp] = -_grad_eta[_qp]/std::sqrt(_norm_sq);
   _dn_dgradphi[_qp] = -I/std::sqrt(_norm_sq) + OG/_norm_cube;
  }
}
