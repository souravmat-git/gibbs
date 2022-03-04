//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeOverall3PStress.h"
registerMooseObject("gibbsApp", ComputeOverall3PStress);

//template <>
InputParameters
ComputeOverall3PStress::validParams()
{
  InputParameters params = Material::validParams();
  return params;
}

ComputeOverall3PStress::ComputeOverall3PStress(const InputParameters & parameters)
  : Material(parameters),
   //Stresses of alpha and beta phase
   _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
   _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
   _gamma_stress(getMaterialProperty<RankTwoTensor>("gamma_stress")),
   //Stiffness of alpha and beta phase
   _alpha_stiffness(getMaterialProperty<RankFourTensor>("alpha_elasticity_tensor")),
   _beta_stiffness(getMaterialProperty<RankFourTensor>("beta_elasticity_tensor")),
   _gamma_stiffness(getMaterialProperty<RankFourTensor>("gamma_elasticity_tensor")),
   //Strain jumps
   _strain_jump_alpha_beta(getMaterialProperty<RankTwoTensor>("strain_jump_alpha_beta")),
   _strain_jump_beta_gamma(getMaterialProperty<RankTwoTensor>("strain_jump_beta_gamma")),
   //and its derivatives wrt to strain and phi
   _ds_alpha_beta_de(getMaterialProperty<RankFourTensor>("ds_alpha_beta_de")),
   _ds_beta_gamma_de(getMaterialProperty<RankFourTensor>("ds_beta_gamma_de")),
   //Derivative of strain jumps with respect to phi_alpha
   _ds_alpha_beta_dphialpha(getMaterialProperty<RankTwoTensor>("ds_alpha_beta_dphialpha")),
   _ds_beta_gamma_dphialpha(getMaterialProperty<RankTwoTensor>("ds_beta_gamma_dphialpha")),
   //Derivative of strain jumps with respect to phi_beta
   _ds_alpha_beta_dphibeta(getMaterialProperty<RankTwoTensor>("ds_alpha_beta_dphibeta")),
   _ds_beta_gamma_dphibeta(getMaterialProperty<RankTwoTensor>("ds_beta_gamma_dphibeta")),
   //Derivative of strain jumps with respect to phi_gamma
   _ds_alpha_beta_dphigamma(getMaterialProperty<RankTwoTensor>("ds_alpha_beta_dphigamma")),
   _ds_beta_gamma_dphigamma(getMaterialProperty<RankTwoTensor>("ds_beta_gamma_dphigamma")),
   //Interpolation functions
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_beta(getMaterialProperty<Real>("h_beta")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //its derivatives with respect to phi_alpha
   _dhalpha_dphialpha(getMaterialProperty<Real>("dhalpha_dphialpha")),
   _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),
   _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
   //its derivatives with respect to phi_beta
   _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
   _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
   //its derivatives with respect to phi_gamma
   _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
   _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma")),
   _dhgamma_dphigamma(getMaterialProperty<Real>("dhgamma_dphigamma")),
   //Compute the following properties
   _sigma_val(declareProperty<RankTwoTensor>("stress")),
   _dsigma_de_val(declareProperty<RankFourTensor>("Jacobian_mult")),
   _dstress_dphialpha(declareProperty<RankTwoTensor>("dstress_dphialpha")),
   _dstress_dphibeta(declareProperty<RankTwoTensor>("dstress_dphibeta")),
   _dstress_dphigamma(declareProperty<RankTwoTensor>("dstress_dphigamma"))
{
}

void
ComputeOverall3PStress::computeQpProperties()
{
    const RankFourTensor avg_C = (_h_alpha[_qp] * _alpha_stiffness[_qp]
                                + _h_beta[_qp]  * _beta_stiffness[_qp]
                                + _h_gamma[_qp] * _gamma_stiffness[_qp]);

    //overall stress
    _sigma_val[_qp] = _alpha_stress[_qp] * _h_alpha[_qp]
                    + _beta_stress[_qp]  * _h_beta[_qp]
                    + _gamma_stress[_qp] * _h_gamma[_qp];

    //Derivative wrt phase alpha variable (Second rank tensor)
    _dstress_dphialpha[_qp] =   _dhbeta_dphialpha[_qp]  * (_beta_stress[_qp]  - _alpha_stress[_qp])
                              + _dhgamma_dphialpha[_qp] * (_gamma_stress[_qp] - _alpha_stress[_qp])
                              + _h_alpha[_qp]* _h_beta[_qp] * (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _ds_alpha_beta_dphialpha[_qp]
                              + _h_beta[_qp] * _h_gamma[_qp] * (_beta_stiffness[_qp] - _gamma_stiffness[_qp])* _ds_beta_gamma_dphialpha[_qp]
                              - _dhalpha_dphialpha[_qp] * avg_C * _strain_jump_alpha_beta[_qp]
                              + _dhgamma_dphialpha[_qp] * avg_C * _strain_jump_beta_gamma[_qp];

    //Derivative wrt phase beta variable (Second rank tensor)
    _dstress_dphibeta[_qp] =   _dhalpha_dphibeta[_qp] * (_alpha_stress[_qp]  - _beta_stress[_qp])
                              + _dhgamma_dphibeta[_qp] * (_gamma_stress[_qp] - _beta_stress[_qp])
                              + _h_alpha[_qp]* _h_beta[_qp] * (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _ds_alpha_beta_dphibeta[_qp]
                              + _h_beta[_qp] * _h_gamma[_qp] * (_beta_stiffness[_qp] - _gamma_stiffness[_qp])* _ds_beta_gamma_dphibeta[_qp]
                              - _dhalpha_dphibeta[_qp] * avg_C * _strain_jump_alpha_beta[_qp]
                              + _dhgamma_dphibeta[_qp] * avg_C * _strain_jump_beta_gamma[_qp];

    //Derivative wrt phase gamma variable (Second rank tensor)
    _dstress_dphigamma[_qp] =   _dhalpha_dphigamma[_qp] * (_alpha_stress[_qp] - _gamma_stress[_qp])
                              + _dhbeta_dphigamma[_qp]  * (_beta_stress[_qp]  - _gamma_stress[_qp])
                              + _h_alpha[_qp]* _h_beta[_qp] * (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _ds_alpha_beta_dphigamma[_qp]
                              + _h_beta[_qp] * _h_gamma[_qp] * (_beta_stiffness[_qp] - _gamma_stiffness[_qp])* _ds_beta_gamma_dphigamma[_qp]
                              - _dhalpha_dphigamma[_qp] * avg_C * _strain_jump_alpha_beta[_qp]
                              + _dhgamma_dphigamma[_qp] * avg_C * _strain_jump_beta_gamma[_qp];

   //Derivative wrt strain
   _dsigma_de_val[_qp] = avg_C
            + _h_alpha[_qp] * _h_beta[_qp] * (_alpha_stiffness[_qp] - _beta_stiffness[_qp]) * _ds_alpha_beta_de[_qp]
            + _h_gamma[_qp] * _h_beta[_qp] * (_beta_stiffness[_qp] - _gamma_stiffness[_qp]) * _ds_beta_gamma_de[_qp];

}
