//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Written by S.Chatterjee

#include "ComputeQDrivingForceOffDiagGamma_3P.h"
registerMooseObject("gibbsApp", ComputeQDrivingForceOffDiagGamma_3P);

template <>
InputParameters
validParams<ComputeQDrivingForceOffDiagGamma_3P>()
{
  InputParameters params = validParams<ComputeQDrivingForceOffDiagBase_3P>();
  return params;
}

ComputeQDrivingForceOffDiagGamma_3P::ComputeQDrivingForceOffDiagGamma_3P(const InputParameters & parameters)
  : ComputeQDrivingForceOffDiagBase_3P(parameters),
   //First derivative of interpolation function
   _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
   _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma")),
   _dhgamma_dphigamma(getMaterialProperty<Real>("dhgamma_dphigamma")),
   //Compute the following properties
   _d2Fdcdstrain_gamma(declareProperty<RankTwoTensor>("d2Fdcdstrain_gamma"))
{
}

void
ComputeQDrivingForceOffDiagGamma_3P::computeQpProperties()
{
  //Driving force with respect to strain for the alpha phase
 _d2Fdcdstrain_gamma[_qp] = _nd_factor[_qp] * (
                        _dhalpha_dphigamma[_qp] * (ComputeQDrivingForceOffDiagBase_3P::dsalpha() - ComputeQDrivingForceOffDiagBase_3P::dsgamma())
                      + _dhbeta_dphigamma[_qp]  * (ComputeQDrivingForceOffDiagBase_3P::dsbeta()  - ComputeQDrivingForceOffDiagBase_3P::dsgamma())
                      - _dhalpha_dphigamma[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_alpha_beta[_qp])
                      - _dhalpha_dphigamma[_qp] * _ds_alpha_beta_de[_qp].innerProductTranspose(_stress[_qp])
                      + _dhgamma_dphigamma[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_beta_gamma[_qp])
                      + _dhgamma_dphigamma[_qp] * _ds_beta_gamma_de[_qp].innerProductTranspose(_stress[_qp]));

}
