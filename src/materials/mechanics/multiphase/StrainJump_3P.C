//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "StrainJump_3P.h"
registerMooseObject("gibbsApp", StrainJump_3P);

//template <>
InputParameters
StrainJump_3P::validParams()
{
  InputParameters params = Material::validParams();
  //params.addRequiredCoupledVar("displacements", "u v w");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_alpha",
  //                "Material constant for alpha phase");
  //params.addRequiredParam<MaterialPropertyName>("stiffness_beta",
  //                "Material constant for alpha phase");
  //params.addRequiredCoupledVar("eta","phase-field variable");
  //params.addRequiredParam<MaterialPropertyName>("total_strain","total_strain_name");
  return params;
}

StrainJump_3P::StrainJump_3P(const InputParameters & parameters)
  : Material(parameters),
  //_gradphi_alpha(coupledGradient("gradphi_alpha")),
  //_gradphi_beta(coupledGradient("gradphi_beta")),
  //_gradphi_gamma(coupledGradient("gradphi_gamm")),
  //Vector array of displacement gradients
  //_grad_disp(coupledGradients("displacements")),
  _total_strain(getMaterialProperty<RankTwoTensor>("total_strain")),
   //Stiffness of alpha and beta phase
   _mat_const_alpha(getMaterialProperty<Real>("alpha_mat_const")),
   _mat_const_beta(getMaterialProperty<Real>("beta_mat_const")),
   _mat_const_gamma(getMaterialProperty<Real>("gamma_mat_const")),
   //Interpolation function
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_beta(getMaterialProperty<Real>("h_beta")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //Compute the following properties
   _a_alpha_beta(declareProperty<Real>("a_alpha_beta")),
   _a_beta_gamma(declareProperty<Real>("a_beta_gamma"))
{
  //unsigned int ndisp = coupledComponents("displacements");
  //Check if number of displacements = mesh dimension
  //if (ndisp != _mesh.dimension())
  //    mooseError(
  //      "The number of variables supplied in 'disp' must match the mesh dimension.");

  //Depending on the mesh dimension the components of disp gradient to zero
  //for(unsigned i = ndisp; i<3; i++)
  //  _grad_disp[i] = &_grad_zero;

}

void
StrainJump_3P::computeQpProperties(){

  //First, calculate the total strain from displacements
  //Second, calculate the difference in elastic constants
  //Third, define the tensors M1, M2, L1, L2
  //Fourth, calculate the inverse tensor P
  //Fifth, calculate a_beta_gamma


  //Define the displacement gradient tensor which is in general not symmetric
  //RankTwoTensor _disp_grad_tensor((*_grad_disp[0])[_qp],
  //                                 (*_grad_disp[1])[_qp],
  //                                 (*_grad_disp[2])[_qp]);


  //For simplicity, we first check the 1D solution
  Real _ex = _total_strain[_qp](0,0);

  //Real _lambda_alpha = _stiffness_alpha[_qp](0,0,1,1);
  //Real _lambda_beta  = _stiffness_beta[_qp](0,0,1,1);
  //Real _lambda_gamma = _stiffness_gamma[_qp](0,0,1,1);

  //mu = C_2323 = C_44 (Voigt Notation)
  //Real  _mu_alpha = _stiffness_alpha[_qp](1,2,1,2);
  //Real  _mu_beta  = _stiffness_alpha[_qp](1,2,1,2);
  //Real  _mu_gamma = _stiffness_gamma[_qp](1,2,1,2);

  //Define mat constant
  //Real _mat_const_alpha = (_lambda_alpha + 2*_mu_alpha);
  //Real _mat_const_beta  = (_lambda_beta  + 2*_mu_beta);
  //Real _mat_const_gamma = (_lambda_gamma + 2*_mu_gamma);

  //Calculate the difference
  Real diff_mc_alphabeta = (_mat_const_alpha[_qp] - _mat_const_beta[_qp]);
  Real diff_mc_betagamma = (_mat_const_beta[_qp]  - _mat_const_gamma[_qp]);

  //Define all prerequiste tensors
  Real L1 = (_h_beta[_qp] + _h_gamma[_qp]) * _mat_const_alpha[_qp] + _h_alpha[_qp] * _mat_const_beta[_qp];
  Real L2 = _h_gamma[_qp] * (_mat_const_alpha[_qp] - _mat_const_beta[_qp]);
  Real M1 = _h_alpha[_qp] * (_mat_const_gamma[_qp] - _mat_const_beta[_qp]);
  Real M2 = _h_gamma[_qp] * _mat_const_beta[_qp] + (_h_beta[_qp] + _h_alpha[_qp])* _mat_const_gamma[_qp];

  Real inv_P = 1.0/(L1*M2 - L2*M1);

  //Calculate the strain jump between alpha/beta
  //and beta/gamma
  _a_alpha_beta[_qp] = -inv_P*(diff_mc_alphabeta * M2 - diff_mc_betagamma * L2)*_ex;
  _a_beta_gamma[_qp] =  inv_P*(diff_mc_alphabeta * M1 - diff_mc_betagamma * L1)*_ex;

}
