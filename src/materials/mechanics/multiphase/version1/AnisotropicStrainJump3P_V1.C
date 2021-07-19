//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AnisotropicStrainJump3P_V1.h"
registerMooseObject("gibbsApp", AnisotropicStrainJump3P_V1);

template <>
InputParameters
validParams<AnisotropicStrainJump3P_V1>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  return params;
}

AnisotropicStrainJump3P_V1::AnisotropicStrainJump3P_V1(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
  //_n1(getMaterialProperty<RealVectorValue>("n1")),
  //_n2(getMaterialProperty<RealVectorValue>("n2")),
  //_total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
   //Stiffness of alpha and beta phase
   _alpha_elasticity_tensor(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _beta_elasticity_tensor(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   _gamma_elasticity_tensor(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
   //Interpolation function
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_beta(getMaterialProperty<Real>("h_beta")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //Obtain the prerequisite tensors
   //_lambda_star(getMaterialProperty<RankTwoTensor>("lambda_star")),
   //_L_hash(getMaterialProperty<RankTwoTensor>("L_hash")),
   _inv_D(getMaterialProperty<RankTwoTensor>("inv_D")),
   _S(getMaterialProperty<RankTwoTensor>("S")),
   //Obtain the prerequisite vectors
   _m_alpha1(getMaterialProperty<RealVectorValue>("m_alpha1")),
   _m_beta1(getMaterialProperty<RealVectorValue>("m_beta1")),
   _psi_beta2(getMaterialProperty<RealVectorValue>("psi_beta2")),
   _psi_gamma2(getMaterialProperty<RealVectorValue>("psi_gamma2")),
   //Compute the following properties
   _a_alpha_beta(declareProperty<RealVectorValue>("avec_alpha_beta")),
   _a_beta_gamma(declareProperty<RealVectorValue>("avec_beta_gamma"))
{
}

void
AnisotropicStrainJump3P_V1::computeQpProperties()
{
  //const Real _un0_norm_sq = _grad_phialpha[_qp].norm_sq();
  //const Real _un1_norm_sq = _grad_phibeta[_qp].norm_sq();
  //const Real _un2_norm_sq = _grad_phigamma[_qp].norm_sq();

  //These are the derived vectors
  RealVectorValue delta_m    = (_m_alpha1[_qp] -  _m_beta1[_qp]);
  RealVectorValue delta_psi  = (_psi_beta2[_qp] - _psi_gamma2[_qp]);

  //Calculate the strain jump between alpha/beta
  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){
      _a_alpha_beta[_qp].zero();
  } else {
      _a_alpha_beta[_qp] = _inv_D[_qp]*(-delta_m);
        //+ _lambda_star[_qp]* _S[_qp] * delta_psi)
   }

   //To calculate beta/gamma we need a_alpha_beta
  if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){
     _a_beta_gamma[_qp].zero();
  } else {
    RealVectorValue b = - (_psi_beta2[_qp] - _psi_gamma2[_qp]);
      //+ _L_hash[_qp]*_a_alpha_beta[_qp]);
    _a_beta_gamma[_qp] = _S[_qp] * b;
  }

}
