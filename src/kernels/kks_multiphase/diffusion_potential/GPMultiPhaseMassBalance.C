//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This was written by S.Chatterjee

#include "GPMultiPhaseMassBalance.h"

registerMooseObject("gibbsApp", GPMultiPhaseMassBalance);

template <>
InputParameters
validParams<GPMultiPhaseMassBalance>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation"
                              "Eqn:  \nabla * (M(\nabla mu) = 0");
  params.addRequiredCoupledVar("phase_alpha", "phase 1");
  params.addRequiredCoupledVar("phase_beta", "Phase 2");
  params.addCoupledVar("phase_gamma", 0.0, "Phase 3");
  params.addRequiredParam<MaterialPropertyName>("h_alpha", "interpolation");
  params.addRequiredParam<MaterialPropertyName>("h_beta", "interpolation");
  params.addParam<MaterialPropertyName>("h_gamma",0.0, "interpolation");
  return params;
}

GPMultiPhaseMassBalance::GPMultiPhaseMassBalance(const InputParameters & parameters)
  :Kernel(parameters),
  _phase_alpha(coupledValue("phase_alpha")),
  _phase_alpha_var(coupled("phase_alpha")),
  _phase_beta(coupledValue("phase_beta")),
  _phase_beta_var(coupled("phase_beta")), 
  //Interpolation Material
  _h_alpha_name(getParam<MaterialPropertyName>("h_alpha")),
  _h_beta_name(getParam<MaterialPropertyName>("h_beta")),
  _h_gamma_name(getParam<MaterialPropertyName>("h_gamma")),
  _h_alpha(getMaterialProperty<Real>(_h_alpha_name)),
  _h_beta(getMaterialProperty<Real>(_h_beta_name)),
  _h_gamma(getMaterialProperty<Real>(_h_gamma_name)),
  // for coupled variable phase_alpha
  _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),    
  // for coupled variable phase_beta
  _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
  //Kinetic material
  _L_BB_beta(getMaterialProperty<Real>("L_BB_beta")),
  _L_BB_alpha(getMaterialProperty<Real>("L_BB_alpha")),
  _dL_BB_muB_beta(getMaterialProperty<Real>("dL_BB_muB_beta")),
  _dL_BB_muB_alpha(getMaterialProperty<Real>("dL_BB_muB_alpha"))
{
}

Real
GPMultiPhaseMassBalance::L_BB_interp() const
{
  return (_L_BB_alpha[_qp]* _h_alpha[_qp] + _L_BB_beta[_qp]*_h_beta[_qp]);
}

Real
GPMultiPhaseMassBalance::dL_BB_muB_interp() const
{
  return (_dL_BB_muB_beta[_qp]* _h_alpha[_qp] + _dL_BB_muB_alpha[_qp]*_h_beta[_qp]);
}

Real
GPMultiPhaseMassBalance::computeQpResidual()
{   
  //Variable on which the kernel operates: mu_B
  return (_grad_test[_i][_qp] * GPMultiPhaseMassBalance::L_BB_interp() * _grad_u[_qp]);
}


Real
GPMultiPhaseMassBalance::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * ( (GPMultiPhaseMassBalance::L_BB_interp() * _grad_phi[_j][_qp])
                               +  (GPMultiPhaseMassBalance::dL_BB_muB_interp() * _grad_u[_qp]) * _phi[_j][_qp]) );
}


Real
GPMultiPhaseMassBalance::computeQpOffDiagJacobian(unsigned int jvar)
{
 if (jvar == _phase_alpha_var)
 {
    return (_grad_test[_i][_qp] *_dhbeta_dphialpha[_qp] * (_L_BB_beta[_qp] - _L_BB_alpha[_qp]) * _grad_u[_qp] * _phi[_j][_qp]);
 }
 else if (jvar == _phase_beta_var)
 {
   return (_grad_test[_i][_qp] *_dhalpha_dphibeta[_qp] * (_L_BB_alpha[_qp] - _L_BB_beta[_qp]) *  _grad_u[_qp] * _phi[_j][_qp]);
 }
  else
  {
    return 0;
  }
}
