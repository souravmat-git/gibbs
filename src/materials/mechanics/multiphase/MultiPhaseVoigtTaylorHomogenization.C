//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*Written by S.Chatterjee*/

#include "MultiPhaseVoigtTaylorHomogenization.h"
registerMooseObject("gibbsApp", MultiPhaseVoigtTaylorHomogenization);

template <>
InputParameters
validParams<MultiPhaseVoigtTaylorHomogenization>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Computes the strain jump [e] tensor for alpha-beta");
  return params;
}

MultiPhaseVoigtTaylorHomogenization::MultiPhaseVoigtTaylorHomogenization(const InputParameters & parameters)
  : Material(parameters),
   //Magnitude of Jump and its derivative wrt h
   _avec_alpha_beta(getMaterialProperty<RealVectorValue>("avec_alpha_beta")),
   _avec_beta_gamma(getMaterialProperty<RealVectorValue>("avec_beta_gamma")),
   //and their derivatives with respect to phi_alpha, phi_beta and phi_gamma
   /*
   _davec_alpha_beta_dphialpha(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphialpha")),
   _davec_alpha_beta_dphibeta(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphibeta")),
   _davec_alpha_beta_dphigamma(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphigamma")),
   //_derivative a_betagamma with respect to phi_alpha, phi_beta and phi_gamma
   _davec_beta_gamma_dphialpha(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphialpha")),
   _davec_beta_gamma_dphibeta(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphibeta")),
   _davec_beta_gamma_dphigamma(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphigamma")),
   */
   //The (compatible) strain jump tensor to be calculated by this class
   _strain_jump_alpha_beta(declareProperty<RankTwoTensor>("strain_jump_alpha_beta")),
   _strain_jump_beta_gamma(declareProperty<RankTwoTensor>("strain_jump_beta_gamma")),
   //Derivative of strain jump with respect to strain to be calculated
   _ds_alpha_beta_de_val(declareProperty<RankFourTensor>("ds_alpha_beta_de")),
   _ds_beta_gamma_de_val(declareProperty<RankFourTensor>("ds_beta_gamma_de")),
   //and its derivatives with respect to phase_alpha
   _ds_alpha_beta_dphialpha(declareProperty<RankTwoTensor>("ds_alpha_beta_dphialpha")),
   _ds_beta_gamma_dphialpha(declareProperty<RankTwoTensor>("ds_beta_gamma_dphialpha")),
   //and its derivatives with respect to phase_beta
   _ds_alpha_beta_dphibeta(declareProperty<RankTwoTensor>("ds_alpha_beta_dphibeta")),
   _ds_beta_gamma_dphibeta(declareProperty<RankTwoTensor>("ds_beta_gamma_dphibeta")),
   //and its derivatives with respect to phase_gamma
   _ds_alpha_beta_dphigamma(declareProperty<RankTwoTensor>("ds_alpha_beta_dphigamma")),
   _ds_beta_gamma_dphigamma(declareProperty<RankTwoTensor>("ds_beta_gamma_dphigamma"))
{
}

void
MultiPhaseVoigtTaylorHomogenization::computeQpProperties(){

     _strain_jump_alpha_beta[_qp].zero();
     _ds_alpha_beta_de_val[_qp].zero();
     _ds_alpha_beta_dphialpha[_qp].zero();
     _ds_alpha_beta_dphibeta[_qp].zero();
     _ds_alpha_beta_dphigamma[_qp].zero();


      _strain_jump_beta_gamma[_qp].zero();
      _ds_beta_gamma_de_val[_qp].zero();
      _ds_beta_gamma_dphialpha[_qp].zero();
      _ds_beta_gamma_dphibeta[_qp].zero();
      _ds_beta_gamma_dphigamma[_qp].zero();


}
