//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AnisotropicStrainJump3P.h"
registerMooseObject("gibbsApp", AnisotropicStrainJump3P);

//template <>
InputParameters
AnisotropicStrainJump3P::validParams()
{
  InputParameters params = VariableUnitNormalBase::validParams();
  return params;
}

AnisotropicStrainJump3P::AnisotropicStrainJump3P(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
   //Obtain the prerequisite tensors
   _inv_D(getMaterialProperty<RankTwoTensor>("inv_D")),
   _S(getMaterialProperty<RankTwoTensor>("S")),
   //Obtain the prerequisite vectors in the alpha-beta region
   _m_alpha1(getMaterialProperty<RealVectorValue>("m_alpha1")),
   _m_beta1(getMaterialProperty<RealVectorValue>("m_beta1")),
   _Z1(getMaterialProperty<RealVectorValue>("Z1")),
   //Obtain the prerequisite vectors in the beta-gamma region
   _psi_beta2(getMaterialProperty<RealVectorValue>("psi_beta2")),
   _psi_gamma2(getMaterialProperty<RealVectorValue>("psi_gamma2")),
   _Z2(getMaterialProperty<RealVectorValue>("Z2")),
   //Compute the following properties
   _a_alpha_beta(declareProperty<RealVectorValue>("avec_alpha_beta")),
   _a_beta_gamma(declareProperty<RealVectorValue>("avec_beta_gamma"))
{
}

void
AnisotropicStrainJump3P::computeQpProperties()
{

  //Calculate the strain jump between alpha/beta
  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){

      _a_alpha_beta[_qp].zero();

  } else {

       RealVectorValue delta_m    = (_m_alpha1[_qp] -  _m_beta1[_qp]);
      _a_alpha_beta[_qp] = _inv_D[_qp]* (_Z1[_qp] - delta_m);

   }

   //To calculate beta/gamma we need a_alpha_beta
  if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0)){

     _a_beta_gamma[_qp].zero();

  } else {

     //We have explicitly assumed that L_hash is zero
     //Since the traction vector at the triple point has no meaning

    RealVectorValue delta_psi = (_psi_beta2[_qp] - _psi_gamma2[_qp]);
    _a_beta_gamma[_qp] = _S[_qp] * (_Z2[_qp] - delta_psi);

  }

}
