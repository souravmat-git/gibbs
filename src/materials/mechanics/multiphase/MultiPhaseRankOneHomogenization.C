//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*Written by S.Chatterjee*/

#include "MultiPhaseRankOneHomogenization.h"
registerMooseObject("gibbsApp", MultiPhaseRankOneHomogenization);

template <>
InputParameters
validParams<MultiPhaseRankOneHomogenization>()
{
  InputParameters params = validParams<VariableUnitNormalBase>();
  params.addClassDescription("Computes the strain jump [e] tensor for alpha-beta");
  return params;
}

MultiPhaseRankOneHomogenization::MultiPhaseRankOneHomogenization(const InputParameters & parameters)
  : VariableUnitNormalBase(parameters),
   //Magnitude of Jump and its derivative wrt h
   _avec_alpha_beta(getMaterialProperty<RealVectorValue>("avec_alpha_beta")),
   _avec_beta_gamma(getMaterialProperty<RealVectorValue>("avec_beta_gamma")),
   //and their derivatives with respect to phi_alpha, phi_beta and phi_gamma
   _davec_alpha_beta_dphialpha(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphialpha")),
   _davec_alpha_beta_dphibeta(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphibeta")),
   _davec_alpha_beta_dphigamma(getMaterialProperty<RealVectorValue>("davec_alpha_beta_dphigamma")),
   //_derivative a_betagamma with respect to phi_alpha, phi_beta and phi_gamma
   _davec_beta_gamma_dphialpha(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphialpha")),
   _davec_beta_gamma_dphibeta(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphibeta")),
   _davec_beta_gamma_dphigamma(getMaterialProperty<RealVectorValue>("davec_beta_gamma_dphigamma")),
   //The (compatible) strain jump tensor to be calculated by this class
   _strain_jump_alpha_beta(declareProperty<RankTwoTensor>("strain_jump_alpha_beta")),
   _strain_jump_beta_gamma(declareProperty<RankTwoTensor>("strain_jump_beta_gamma")),
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
MultiPhaseRankOneHomogenization::computeQpProperties(){

  if (unalpha_norm_sq() < std::pow(libMesh::TOLERANCE,3.0) ||
      unbeta_norm_sq()  < std::pow(libMesh::TOLERANCE,3.0)){

     _strain_jump_alpha_beta[_qp].zero();
     _ds_alpha_beta_dphialpha[_qp].zero();
     _ds_alpha_beta_dphibeta[_qp].zero();
     _ds_alpha_beta_dphigamma[_qp].zero();

 } else {

   RealVectorValue _n1 = -_grad_phibeta[_qp]/std::sqrt(unbeta_norm_sq());

   RankTwoTensor _jump_in_disp_grad_alpha_beta,
                 _dnablau_alpha_beta_dphialpha,
                 _dnablau_alpha_beta_dphibeta,
                 _dnablau_alpha_beta_dphigamma;

    //First take the vector outer product
   _jump_in_disp_grad_alpha_beta.vectorOuterProduct(_avec_alpha_beta[_qp], _n1);
   _dnablau_alpha_beta_dphialpha.vectorOuterProduct(_davec_alpha_beta_dphialpha[_qp], _n1);
   _dnablau_alpha_beta_dphibeta.vectorOuterProduct(_davec_alpha_beta_dphibeta[_qp], _n1);
   _dnablau_alpha_beta_dphigamma.vectorOuterProduct(_davec_alpha_beta_dphigamma[_qp], _n1);

    //Then use the definition of strain tensor
   _strain_jump_alpha_beta[_qp] = (_jump_in_disp_grad_alpha_beta
                     +  _jump_in_disp_grad_alpha_beta.transpose()) / 2.0;

   _ds_alpha_beta_dphialpha[_qp] = (_dnablau_alpha_beta_dphialpha
                     +  _dnablau_alpha_beta_dphialpha.transpose()) / 2.0;

   _ds_alpha_beta_dphibeta[_qp] = (_dnablau_alpha_beta_dphibeta
                     +  _dnablau_alpha_beta_dphibeta.transpose()) / 2.0;

   _ds_alpha_beta_dphigamma[_qp] = (_dnablau_alpha_beta_dphigamma
                     +  _dnablau_alpha_beta_dphigamma.transpose()) / 2.0;

 }


  //Calculation is restricted to beta-gamma interface
  if (ungamma_norm_sq() < std::pow(libMesh::TOLERANCE,3.0))
  {
      _strain_jump_beta_gamma[_qp].zero();
      _ds_beta_gamma_dphialpha[_qp].zero();
      _ds_beta_gamma_dphibeta[_qp].zero();
      _ds_beta_gamma_dphigamma[_qp].zero();

  } else {

   RealVectorValue _n2 = -_grad_phigamma[_qp]/std::sqrt(ungamma_norm_sq());

   RankTwoTensor _jump_in_disp_grad_beta_gamma,
                 _dnablau_beta_gamma_dphialpha,
                 _dnablau_beta_gamma_dphibeta,
                 _dnablau_beta_gamma_dphigamma;

   //n2 is the unit normal for beta-gamma interface
   _jump_in_disp_grad_beta_gamma.vectorOuterProduct(_avec_beta_gamma[_qp], _n2);
   _dnablau_beta_gamma_dphialpha.vectorOuterProduct(_davec_beta_gamma_dphialpha[_qp], _n2);
   _dnablau_beta_gamma_dphibeta.vectorOuterProduct(_davec_beta_gamma_dphibeta[_qp], _n2);
   _dnablau_beta_gamma_dphigamma.vectorOuterProduct(_davec_beta_gamma_dphigamma[_qp], _n2);

   _strain_jump_beta_gamma[_qp] = (_jump_in_disp_grad_beta_gamma
                                +  _jump_in_disp_grad_beta_gamma.transpose()) / 2.0;

   _ds_beta_gamma_dphialpha[_qp] = (_dnablau_beta_gamma_dphialpha
                                 +  _dnablau_beta_gamma_dphialpha.transpose()) / 2.0;

    _ds_beta_gamma_dphibeta[_qp] = (_dnablau_beta_gamma_dphibeta
                                 +  _dnablau_beta_gamma_dphibeta.transpose()) / 2.0;

    _ds_beta_gamma_dphigamma[_qp] = (_dnablau_beta_gamma_dphigamma
                                  +  _dnablau_beta_gamma_dphigamma.transpose()) / 2.0;

  }

}
