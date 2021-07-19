//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Written by S.Chatterjee

#include "ComputeQDrivingForceOffDiagAlpha_3P.h"
registerMooseObject("gibbsApp", ComputeQDrivingForceOffDiagAlpha_3P);

template <>
InputParameters
validParams<ComputeQDrivingForceOffDiagAlpha_3P>()
{
  InputParameters params = validParams<ComputeQDrivingForceOffDiagBase_3P>();
  return params;
}

ComputeQDrivingForceOffDiagAlpha_3P::ComputeQDrivingForceOffDiagAlpha_3P(const InputParameters & parameters)
  : ComputeQDrivingForceOffDiagBase_3P(parameters),
   //First derivative of interpolation function
   _dhalpha_dphialpha(getMaterialProperty<Real>("dhalpha_dphialpha")),
   _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),
   _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
   //Compute the following properties
   _d2Fdcdstrain_alpha(declareProperty<RankTwoTensor>("d2Fdcdstrain_alpha"))
{
}

void
ComputeQDrivingForceOffDiagAlpha_3P::computeQpProperties()
{
  //Driving force with respect to strain for the alpha phase
 _d2Fdcdstrain_alpha[_qp] = _nd_factor[_qp] * (
                  _dhbeta_dphialpha[_qp]  * (ComputeQDrivingForceOffDiagBase_3P::dsbeta() - ComputeQDrivingForceOffDiagBase_3P::dsalpha())
                + _dhgamma_dphialpha[_qp] * (ComputeQDrivingForceOffDiagBase_3P::dsgamma() - ComputeQDrivingForceOffDiagBase_3P::dsalpha())
                - _dhalpha_dphialpha[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_alpha_beta[_qp])
                - _dhalpha_dphialpha[_qp] * _ds_alpha_beta_de[_qp].innerProductTranspose(_stress[_qp])
                + _dhgamma_dphialpha[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_beta_gamma[_qp])
                + _dhgamma_dphialpha[_qp] * _ds_beta_gamma_de[_qp].innerProductTranspose(_stress[_qp]));

}
