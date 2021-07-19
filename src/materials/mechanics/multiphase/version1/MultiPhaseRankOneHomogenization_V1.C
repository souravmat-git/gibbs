//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

/*Written by S.Chatterjee*/

#include "MultiPhaseRankOneHomogenization_V1.h"
registerMooseObject("gibbsApp", MultiPhaseRankOneHomogenization_V1);

template <>
InputParameters
validParams<MultiPhaseRankOneHomogenization_V1>()
{
  InputParameters params = validParams<Material>();
  //params.addRequiredParam<MaterialPropertyName>("strain_jump_name",
  //                          "Name of strain jump for alpha");
  params.addRequiredCoupledVar("phase_alpha", "alpha phase");
  params.addRequiredCoupledVar("phase_beta",  "beta phase");
  params.addRequiredCoupledVar("phase_gamma", "gamma phase");
  params.addClassDescription("Computes the strain jump [e] tensor for alpha-beta");
  return params;
}

MultiPhaseRankOneHomogenization_V1::MultiPhaseRankOneHomogenization_V1(const InputParameters & parameters)
  : Material(parameters),
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
   //Fetch the variables
   _grad_alpha(coupledGradient("phase_alpha")),
   _grad_beta(coupledGradient("phase_beta")),
   _grad_gamma(coupledGradient("phase_gamma")),
   //Unit normals
   _n1(getMaterialProperty<RealVectorValue>("n1")),
   _n2(getMaterialProperty<RealVectorValue>("n2")),
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
MultiPhaseRankOneHomogenization_V1::computeQpProperties(){

  //Calculation is restricted to alpha-beta interface
  //if (_norm_sq_alpha < tol && _norm_sq_beta < tol){
  //   _strain_jump_alpha_beta[_qp].zero();
  //}
  //else {

  RankTwoTensor _jump_in_disp_grad_alpha_beta,
                _dnablau_alpha_beta_dphialpha,
                _dnablau_alpha_beta_dphibeta,
                _dnablau_alpha_beta_dphigamma;

  _jump_in_disp_grad_alpha_beta.vectorOuterProduct(_avec_alpha_beta[_qp], _n1[_qp]);
  _dnablau_alpha_beta_dphialpha.vectorOuterProduct(_davec_alpha_beta_dphialpha[_qp], _n1[_qp]);
  _dnablau_alpha_beta_dphibeta.vectorOuterProduct(_davec_alpha_beta_dphibeta[_qp], _n1[_qp]);
  _dnablau_alpha_beta_dphigamma.vectorOuterProduct(_davec_alpha_beta_dphigamma[_qp], _n1[_qp]);

  _strain_jump_alpha_beta[_qp] = (_jump_in_disp_grad_alpha_beta
                    +  _jump_in_disp_grad_alpha_beta.transpose()) / 2.0;

  _ds_alpha_beta_dphialpha[_qp] = (_dnablau_alpha_beta_dphialpha
                    +  _dnablau_alpha_beta_dphialpha.transpose()) / 2.0;

  _ds_alpha_beta_dphibeta[_qp] = (_dnablau_alpha_beta_dphibeta
                    +  _dnablau_alpha_beta_dphibeta.transpose()) / 2.0;

  _ds_alpha_beta_dphigamma[_qp] = (_dnablau_alpha_beta_dphigamma
                    +  _dnablau_alpha_beta_dphigamma.transpose()) / 2.0;


  //Calculation is restricted to beta-gamma interface
  //if (_norm_sq_beta < tol && _norm_sq_gamma < tol){
  //   _strain_jump_beta_gamma[_qp].zero();
  //}
  //else{

   RankTwoTensor _jump_in_disp_grad_beta_gamma,
                 _dnablau_beta_gamma_dphialpha,
                 _dnablau_beta_gamma_dphibeta,
                 _dnablau_beta_gamma_dphigamma;

   //n2 is the unit normal for beta-gamma interface
  _jump_in_disp_grad_beta_gamma.vectorOuterProduct(_avec_beta_gamma[_qp], _n2[_qp]);
  _dnablau_beta_gamma_dphialpha.vectorOuterProduct(_davec_beta_gamma_dphialpha[_qp], _n2[_qp]);
  _dnablau_beta_gamma_dphibeta.vectorOuterProduct(_davec_beta_gamma_dphibeta[_qp], _n2[_qp]);
  _dnablau_beta_gamma_dphigamma.vectorOuterProduct(_davec_beta_gamma_dphigamma[_qp], _n2[_qp]);

  _strain_jump_beta_gamma[_qp] = (_jump_in_disp_grad_beta_gamma
                               +  _jump_in_disp_grad_beta_gamma.transpose()) / 2.0;

  _ds_beta_gamma_dphialpha[_qp] = (_dnablau_beta_gamma_dphialpha
                                + _dnablau_beta_gamma_dphialpha.transpose()) / 2.0;

  _ds_beta_gamma_dphibeta[_qp] = (_dnablau_beta_gamma_dphibeta
                               +  _dnablau_beta_gamma_dphibeta.transpose()) / 2.0;

  _ds_beta_gamma_dphigamma[_qp] = (_dnablau_beta_gamma_dphigamma
                                +  _dnablau_beta_gamma_dphigamma.transpose()) / 2.0;

}
