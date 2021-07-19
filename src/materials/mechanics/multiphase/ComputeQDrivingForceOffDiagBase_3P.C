//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//

#include "ComputeQDrivingForceOffDiagBase_3P.h"
registerMooseObject("gibbsApp", ComputeQDrivingForceOffDiagBase_3P);

template <>
InputParameters
validParams<ComputeQDrivingForceOffDiagBase_3P>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<MaterialPropertyName>("nd_factor", "Non-dimensional factor");
  return params;
}

ComputeQDrivingForceOffDiagBase_3P::ComputeQDrivingForceOffDiagBase_3P(const InputParameters & parameters)
  : Material(parameters),
   //Stresses of alpha, beta and gamma phases
   _alpha_stress(getMaterialProperty<RankTwoTensor>("alpha_stress")),
   _beta_stress(getMaterialProperty<RankTwoTensor>("beta_stress")),
   _gamma_stress(getMaterialProperty<RankTwoTensor>("gamma_stress")),
   //Strain jumps
   _strain_jump_alpha_beta(getMaterialProperty<RankTwoTensor>("strain_jump_alpha_beta")),
   _strain_jump_beta_gamma(getMaterialProperty<RankTwoTensor>("strain_jump_beta_gamma")),
   //interpolation function
   _h_alpha(getMaterialProperty<Real>("h_alpha")),
   _h_beta(getMaterialProperty<Real>("h_beta")),
   _h_gamma(getMaterialProperty<Real>("h_gamma")),
   //Derivative of strain jumps with respect to strain
   _ds_alpha_beta_de(getMaterialProperty<RankFourTensor>("ds_alpha_beta_de")),
   _ds_beta_gamma_de(getMaterialProperty<RankFourTensor>("ds_beta_gamma_de")),
   //overall stress
   _stress(getMaterialProperty<RankTwoTensor>("stress")),
   //Derivative of overall stress with respect to strain
   _Jacobian_mult(getMaterialProperty<RankFourTensor>("Jacobian_mult")),
   _nd_factor(getMaterialProperty<Real>("nd_factor"))
{
}

//We first calculate these tensors required in derived classes
RankTwoTensor
ComputeQDrivingForceOffDiagBase_3P::dsalpha() const{

  return  _alpha_stress[_qp]
         + (_h_beta[_qp] + _h_gamma[_qp])* _ds_alpha_beta_de[_qp].innerProductTranspose(_alpha_stress[_qp])
         + _h_gamma[_qp]* _ds_beta_gamma_de[_qp].innerProductTranspose(_alpha_stress[_qp]);
}

RankTwoTensor
ComputeQDrivingForceOffDiagBase_3P::dsbeta() const{

  return _beta_stress[_qp]
        - _h_alpha[_qp] * _ds_alpha_beta_de[_qp].innerProductTranspose(_beta_stress[_qp])
        + _h_gamma[_qp] * _ds_beta_gamma_de[_qp].innerProductTranspose(_beta_stress[_qp]);
}

RankTwoTensor
ComputeQDrivingForceOffDiagBase_3P::dsgamma() const{

   return _gamma_stress[_qp]
          - _h_alpha[_qp] * _ds_alpha_beta_de[_qp].innerProductTranspose(_gamma_stress[_qp])
          - (_h_beta[_qp] + _h_alpha[_qp]) * _ds_beta_gamma_de[_qp].innerProductTranspose(_gamma_stress[_qp]);
}


void
ComputeQDrivingForceOffDiagBase_3P::computeQpProperties(){}
