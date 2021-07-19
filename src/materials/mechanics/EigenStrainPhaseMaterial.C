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

template <>
InputParameters
validParams<EigenStrainPhaseMaterial>()
{
  InputParameters params = validParams<ComputeEigenstrain>();
  return params;
}

EigenStrainPhaseMaterial::EigenStrainPhaseMaterial(const InputParameters & parameters)
  : ComputeEigenstrain(parameters),
  //Note that in ComputeEigenStrain we obtain _base_name + _
   _dprefactor(getMaterialProperty<Real>(_base_name + "dprefactor")),
   _d2prefactor(getMaterialProperty<Real>(_base_name + "d2prefactor")),
   _deigen_dmu(declareProperty<RankTwoTensor>(_base_name + "deigen_dmu")),
   _d2eigen_dmu2(declareProperty<RankTwoTensor>(_base_name + "d2eigen_dmu2"))
{
} 

void
EigenStrainPhaseMaterial::computeQpEigenstrain(){

    ComputeEigenstrain::computeQpEigenstrain();
    
   _deigen_dmu[_qp]   =  _eigen_base_tensor * _dprefactor[_qp];
   _d2eigen_dmu2[_qp] =  _eigen_base_tensor * _d2prefactor[_qp];

}
