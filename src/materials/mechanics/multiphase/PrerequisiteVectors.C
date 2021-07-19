//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PrerequisiteVectors.h"
registerMooseObject("gibbsApp", PrerequisiteVectors);

template <>
InputParameters
validParams<PrerequisiteVectors>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  //params.addRequiredCoupledVar("displacements", "u v w");
  params.addClassDescription("This material calculates the prerequsite"
                              "second-rank tensors for multiphase case");
  return params;
}

PrerequisiteVectors::PrerequisiteVectors(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
  //total strain
  _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
  //fetch the eigen strains
  _alpha_eigen(getMaterialProperty<RankTwoTensor>("alpha_eigen")),
  _beta_eigen(getMaterialProperty<RankTwoTensor>("beta_eigen")),
  _gamma_eigen(getMaterialProperty<RankTwoTensor>("gamma_eigen")),
  //Stiffness of alpha, beta and gamma phases
  _alpha_elasticity_tensor(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
  _beta_elasticity_tensor(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
  _gamma_elasticity_tensor(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
  //Compute the following properties in the alpha-beta region
  _m_alpha1(declareProperty<RealVectorValue>("m_alpha1")),
  _m_beta1(declareProperty<RealVectorValue>("m_beta1")),
  _Z1(declareProperty<RealVectorValue>("Z1")),
  //Compute the following properties in the beta-gamma regions
  _psi_beta2(declareProperty<RealVectorValue>("psi_beta2")),
  _psi_gamma2(declareProperty<RealVectorValue>("psi_gamma2")),
  _Z2(declareProperty<RealVectorValue>("Z2"))
{
}

void
PrerequisiteVectors::computeQpProperties(){

  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){

      _m_alpha1[_qp].zero();
      _m_beta1[_qp].zero();
      _Z1[_qp].zero();

  } else {

     //Define the unit normal at the alpha-beta interface
     RealVectorValue _n1 = - _grad_phibeta[_qp]/std::sqrt(unbeta_norm_sq());

     //Run the calculation for vectors m_alpha1 and m_beta1
    _m_alpha1[_qp]   = _alpha_elasticity_tensor[_qp] * _total_strain[_qp] * _n1;
    _m_beta1[_qp]    = _beta_elasticity_tensor[_qp]  * _total_strain[_qp] * _n1;
    _Z1[_qp]         = (_alpha_elasticity_tensor[_qp] * _alpha_eigen[_qp]
                     -  _beta_elasticity_tensor[_qp]  * _beta_eigen[_qp])* _n1;

  }

  if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){

      _psi_beta2[_qp].zero();
      _psi_gamma2[_qp].zero();
      _Z2[_qp].zero();

  } else {

    //Define the unit normal at beta-gamma interface
    RealVectorValue _n2 = -_grad_phigamma[_qp]/std::sqrt(ungamma_norm_sq());

    _psi_beta2[_qp]  = _beta_elasticity_tensor[_qp]  * _total_strain[_qp] * _n2;
    _psi_gamma2[_qp] = _gamma_elasticity_tensor[_qp] * _total_strain[_qp] * _n2;
    _Z2[_qp]         = (_beta_elasticity_tensor[_qp] * _beta_eigen[_qp]
                     -  _gamma_elasticity_tensor[_qp] * _gamma_eigen[_qp]) * _n2;

  }

}
