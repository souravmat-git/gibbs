//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This was written by S.Chatterjee

#include "MultiCompContinuityEquationC.h"

registerMooseObject("gibbsApp", MultiCompContinuityEquationC);

template <>
InputParameters
validParams<MultiCompContinuityEquationC>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("This kernel implements the continuity equation on mu_{1}"
                              "Eqn:  nabla * (M_{11}(nabla mu_{1}) + nabla * (M_{12}(nabla mu_{2})) = 0");
  params.addRequiredCoupledVar("B_diff_pot", "diffusion potential of component B");
  params.addRequiredCoupledVar("eta", "Phase field variable");
  params.addCoupledVar("xB_beta", 0.0, "Component B in beta phase");
  params.addCoupledVar("xC_beta", 0.0, "Component C in beta phase");
  params.addCoupledVar("xB_alpha", 0.0, "Component B in alpha phase");
  params.addCoupledVar("xC_alpha", 0.0, "Component C in alpha phase");
  params.addCoupledVar("D_diff_pot", 0.0,"Diffusion potential of comp D");
  params.addParam<MaterialPropertyName>("L_CD_beta",0.0,"Effect of muD on xC");
  params.addParam<MaterialPropertyName>("L_CD_alpha", 0.0, "Effect of mu(D) on xC");
   params.addParam<MaterialPropertyName>("dL_CC_xB_beta", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_BC_xB_beta", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_CC_xC_beta", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_BC_xC_beta", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_CC_xB_alpha", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_BC_xB_alpha", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_CC_xC_alpha", 0.0, "");
  params.addParam<MaterialPropertyName>("dL_BC_xC_alpha", 0.0, "");
  return params;
}

MultiCompContinuityEquationC::MultiCompContinuityEquationC(const InputParameters & parameters)
  : Kernel(parameters),
    _grad_B_diff_pot(coupledGradient("B_diff_pot")),
    _grad_B_diff_pot_var(coupled("B_diff_pot")),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    //Component B & C in beta phase
    _xB_beta(coupledValue("xB_beta")),
    _xB_beta_var(coupled("xB_beta")),
    _xC_beta(coupledValue("xC_beta")),
    _xC_beta_var(coupled("xC_beta")),
    //Component B & C in alpha phase
    _xB_alpha(coupledValue("xB_alpha")),
    _xB_alpha_var(coupled("xB_alpha")),
    _xC_alpha(coupledValue("xC_alpha")),
    _xC_alpha_var(coupled("xC_alpha")),
    _grad_D_diff_pot(coupledGradient("D_diff_pot")),
    _grad_D_diff_pot_var(coupled("D_diff_pot")),
    _L_CC_beta(getMaterialProperty<Real>("L_CC_beta")),
    _L_BC_beta(getMaterialProperty<Real>("L_BC_beta")),
    _L_CD_beta(getMaterialProperty<Real>("L_CD_beta")),
    _L_CC_alpha(getMaterialProperty<Real>("L_CC_alpha")),
    _L_BC_alpha(getMaterialProperty<Real>("L_BC_alpha")),
    _L_CD_alpha(getMaterialProperty<Real>("L_CD_alpha")),
     //Derivatives of L with respect to xB_beta
    _dL_CC_xB_beta(getMaterialProperty<Real>("dL_CC_xB_beta")),
    _dL_BC_xB_beta(getMaterialProperty<Real>("dL_BC_xB_beta")),
    //Derivative of L with respect to xC_beta
    _dL_CC_xC_beta(getMaterialProperty<Real>("dL_CC_xC_beta")),
    _dL_BC_xC_beta(getMaterialProperty<Real>("dL_BC_xC_beta")),
    //Derivative of L with respect to xB_alpha
    _dL_CC_xB_alpha(getMaterialProperty<Real>("dL_CC_xB_alpha")),
    _dL_BC_xB_alpha(getMaterialProperty<Real>("dL_BC_xB_alpha")),
    //Derivative of L with respect to xC_alpha
    _dL_CC_xC_alpha(getMaterialProperty<Real>("dL_CC_xC_alpha")),
    _dL_BC_xC_alpha(getMaterialProperty<Real>("dL_BC_xC_alpha")),    
    _h(getMaterialProperty<Real>("h")),
    _dh(getMaterialProperty<Real>("dh"))
{
}

Real
MultiCompContinuityEquationC::computeQpResidual()
{
  // $nabla M_{11} * nabla mu_{1} + nabla M_{12} nabla mu_{2}$
  
  const Real _L_CC_interp = _L_CC_beta[_qp]* _h[_qp] + _L_CC_alpha[_qp]*(1.0 - _h[_qp]);
  const Real _L_BC_interp = _L_BC_beta[_qp]* _h[_qp] + _L_BC_alpha[_qp]*(1.0 - _h[_qp]);
  const Real _L_CD_interp = _L_CD_beta[_qp]* _h[_qp] + _L_CD_alpha[_qp]*(1.0 - _h[_qp]);


  return _grad_test[_i][_qp] * (_L_CC_interp * _grad_u[_qp] 
                              + _L_BC_interp * _grad_B_diff_pot[_qp]
                              + _L_CD_interp * _grad_D_diff_pot[_qp]);
}

Real
MultiCompContinuityEquationC::computeQpJacobian()
{
  return (_grad_test[_i][_qp] * (_L_CC_beta[_qp]* _h[_qp]
                              +  _L_CC_alpha[_qp] * (1.0 -_h[_qp])) * _grad_phi[_j][_qp]);
}

Real
MultiCompContinuityEquationC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _grad_B_diff_pot_var)
  {
    return (_grad_test[_i][_qp] *(_L_BC_beta[_qp] * _h[_qp]
                                + _L_BC_alpha[_qp] * (1.0 - _h[_qp]))  * _grad_phi[_j][_qp]);
  }
  else if (jvar == _grad_D_diff_pot_var)
  {
    return (_grad_test[_i][_qp] *(_L_CD_beta[_qp] * _h[_qp]
                               +  _L_CD_alpha[_qp] * (1.0 - _h[_qp])) * _grad_phi[_j][_qp]);
  }
  else if (jvar == _eta_var)
  {
    return (_grad_test[_i][_qp] * _dh[_qp] * ((_L_CC_beta[_qp] - _L_CC_alpha[_qp]) * _grad_u[_qp]
                                            + (_L_BC_beta[_qp] - _L_BC_alpha[_qp]) * _grad_B_diff_pot[_qp]
                                            + (_L_CD_beta[_qp] - _L_CD_alpha[_qp]) * _grad_D_diff_pot[_qp])*_phi[_j][_qp]);
  }
    else if (jvar == _xB_beta_var)
  {
    return (_grad_test[_i][_qp] * _h[_qp] * 
                (_dL_CC_xB_beta[_qp] * _grad_u[_qp] + _dL_BC_xB_beta[_qp] * _grad_B_diff_pot[_qp]) *_phi[_j][_qp]); 
  }
  else if (jvar == _xC_beta_var)
  {
    return (_grad_test[_i][_qp] * _h[_qp] *
                (_dL_CC_xC_beta[_qp] * _grad_u[_qp] + _dL_BC_xC_beta[_qp] * _grad_B_diff_pot[_qp]) * _phi[_j][_qp]);
  }
  else if (jvar == _xB_alpha_var)
  {
    return (_grad_test[_i][_qp] * _h[_qp] *
                 (_dL_CC_xB_alpha[_qp] * _grad_u[_qp] + _dL_BC_xB_alpha[_qp] * _grad_B_diff_pot[_qp]) * _phi[_j][_qp]);
  }
  else if (jvar == _xC_alpha_var)
  {
    return (_grad_test[_i][_qp] * _h[_qp] *
                  (_dL_CC_xC_alpha[_qp] * _grad_u[_qp] + _dL_BC_xC_alpha[_qp] * _grad_B_diff_pot[_qp]) * _phi[_j][_qp]);
  }
  else
  {
    return 0.0; 
  }
}
