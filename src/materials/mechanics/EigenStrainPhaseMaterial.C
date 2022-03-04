//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//* Calculates the elastic constants for a given youngs modulus and poissons

#include "EigenStrainPhaseMaterial.h"
registerMooseObject("gibbsApp", EigenStrainPhaseMaterial);

//template <>
InputParameters
EigenStrainPhaseMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<std::string>("base_name", "Phase name");
  //params.addRequiredParam<MaterialPropertyName>("eigen_name",  "eigen")
  return params;
}

EigenStrainPhaseMaterial::EigenStrainPhaseMaterial(const InputParameters & parameters)
  : Material(parameters),
   _base_name(getParam<std::string>("base_name")),
  //Note that in ComputeEigenStrain we obtain _base_name + _
   _eigen_tensor(getMaterialProperty<RankTwoTensor>(_base_name + "_eigen")),
   _dprefactor(getMaterialProperty<Real>(_base_name + "_dprefactor")),
   _d2prefactor(getMaterialProperty<Real>(_base_name + "_d2prefactor")),
   _deigen_dmu(declareProperty<RankTwoTensor>(_base_name + "_deigen_dmu")),
   _d2eigen_dmu2(declareProperty<RankTwoTensor>(_base_name + "_d2eigen_dmu2"))
{
}

void
EigenStrainPhaseMaterial::computeQpProperties(){

    //ComputeEigenstrain::computeQpEigenstrain();

   _deigen_dmu[_qp]   =  _eigen_tensor[_qp] * _dprefactor[_qp];
   _d2eigen_dmu2[_qp] =  _eigen_tensor[_qp] * _d2prefactor[_qp];
}
