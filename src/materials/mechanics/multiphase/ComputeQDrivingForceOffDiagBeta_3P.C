//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Written by S.Chatterjee

#include "ComputeQDrivingForceOffDiagBeta_3P.h"
registerMooseObject("gibbsApp", ComputeQDrivingForceOffDiagBeta_3P);

template <>
InputParameters
validParams<ComputeQDrivingForceOffDiagBeta_3P>()
{
  InputParameters params = validParams<ComputeQDrivingForceOffDiagBase_3P>();
  return params;
}

ComputeQDrivingForceOffDiagBeta_3P::ComputeQDrivingForceOffDiagBeta_3P(const InputParameters & parameters)
  : ComputeQDrivingForceOffDiagBase_3P(parameters),
   //First derivative of interpolation function
   _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
   _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
   //Compute the following properties
   _d2Fdcdstrain_beta(declareProperty<RankTwoTensor>("d2Fdcdstrain_beta"))
{
}

void
ComputeQDrivingForceOffDiagBeta_3P::computeQpProperties()
{
  //Driving force with respect to strain for the alpha phase
 _d2Fdcdstrain_beta[_qp] = _nd_factor[_qp] * (
                        _dhalpha_dphibeta[_qp]  * (ComputeQDrivingForceOffDiagBase_3P::dsalpha() - ComputeQDrivingForceOffDiagBase_3P::dsbeta())
                      + _dhgamma_dphibeta[_qp]  * (ComputeQDrivingForceOffDiagBase_3P::dsgamma() - ComputeQDrivingForceOffDiagBase_3P::dsbeta())
                      - _dhalpha_dphibeta[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_alpha_beta[_qp])
                      - _dhalpha_dphibeta[_qp] * _ds_alpha_beta_de[_qp].innerProductTranspose(_stress[_qp])
                      + _dhgamma_dphibeta[_qp] * _Jacobian_mult[_qp].innerProductTranspose(_strain_jump_beta_gamma[_qp])
                      + _dhgamma_dphibeta[_qp] * _ds_beta_gamma_de[_qp].innerProductTranspose(_stress[_qp]));

}
