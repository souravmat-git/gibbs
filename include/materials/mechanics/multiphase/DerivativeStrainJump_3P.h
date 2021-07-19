//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class DerivativeStrainJump_3P;

//MOOSe includes
#include "VariableUnitNormalBase.h"

template <>
InputParameters validParams<DerivativeStrainJump_3P>();

class DerivativeStrainJump_3P : public VariableUnitNormalBase
{
public:
  DerivativeStrainJump_3P(const InputParameters & parameters);

   //Calculates the derivative of strain jump with respect to strain
   //and phase-field variable

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

  //Strain jumps
  const MaterialProperty<RealVectorValue> & _avec_alpha_beta;
  const MaterialProperty<RealVectorValue> & _psi_beta2;
  const MaterialProperty<RealVectorValue> & _psi_gamma2;

  //Fetch the prerequiste second rank tensors
  const MaterialProperty<RankTwoTensor> & _inv_D;
  const MaterialProperty<RankTwoTensor> & _S;

  //Fetch the stiffness tensors
  const MaterialProperty<RankFourTensor> & _alpha_stiffness;
  const MaterialProperty<RankFourTensor> & _beta_stiffness;
  const MaterialProperty<RankFourTensor> & _gamma_stiffness;

  //Derivative of D with respect to phi_alpha, phi_beta, phi_gamma
  const MaterialProperty<RankTwoTensor>  & _dD_dphialpha;
  const MaterialProperty<RankTwoTensor>  & _dD_dphibeta;
  const MaterialProperty<RankTwoTensor>  & _dD_dphigamma;

  //Derivative of S with resect to phi_alpha, phi_beta and phi_gamma
  const MaterialProperty<RankTwoTensor>  & _dS_dphialpha;
  const MaterialProperty<RankTwoTensor>  & _dS_dphibeta;
  const MaterialProperty<RankTwoTensor>  & _dS_dphigamma;

  //Derivatives of strain jumps with resect to strain
  MaterialProperty<RankFourTensor> & _ds_alpha_beta_de_val;
  MaterialProperty<RankFourTensor> & _ds_beta_gamma_de_val;

  //Derivative of a_alpha_beta with phi_alpha, phi_beta and phi_gamma
  MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphialpha_val;
  MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphibeta_val;
  MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphigamma_val;

  //Derivative of a_betagamma with phi_alpha, phi_beta and phi_gamma
  MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphialpha_val;
  MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphibeta_val;
  MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphigamma_val;
};
