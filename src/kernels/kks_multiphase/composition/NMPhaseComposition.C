//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This code was written by S.Chatterjee

#include "NMPhaseComposition.h"
registerMooseObject("gibbsApp", NMPhaseComposition);

//*class template specialization
template <>
InputParameters
validParams<NMPhaseComposition>()
{
    
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: c = c_{alpha}h_{alpha} + ...\
                              +c_{beta}h_{beta} + c_{gamma}h_{gamma}");
  //The non-variable that this kernel operates on is c_{gamma} 
  params.addCoupledVar("phase_comp_beta", 0.0,"Phase comp. in the beta phase");
  params.addCoupledVar("phase_comp_alpha",0.0,"Phase comp. in the alpha phase");
  params.addCoupledVar("phase_alpha", "Phase field representing the alpha phase");
  params.addCoupledVar("phase_beta",0.0,"Phase field for beta phase");
  params.addCoupledVar("phase_gamma",0.0, "Phase field for the gamma phase");                          
  params.addRequiredCoupledVar("mole_fraction","mole fraction");      
  return params;
}

NMPhaseComposition::NMPhaseComposition(const InputParameters & parameters)
  : Kernel(parameters),
    _phase_comp_beta(coupledValue("phase_comp_beta")),
    _phase_comp_beta_var(coupled("phase_comp_beta")),
    _phase_comp_alpha(coupledValue("phase_comp_alpha")),
    _phase_comp_alpha_var(coupled("phase_comp_alpha")),
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_alpha_var(coupled("phase_alpha")),
    _phase_beta(coupledValue("phase_beta")),
    _phase_beta_var(coupled("phase_beta")),
    _phase_gamma(coupledValue("phase_gamma")),
    _phase_gamma_var(coupled("phase_gamma")),
    _mole_fraction(coupledValue("mole_fraction")),
    _mole_fraction_var(coupled("mole_fraction")),
    _h_alpha(getMaterialProperty<Real>("h_alpha")),
    _h_beta(getMaterialProperty<Real>("h_beta")),
    _h_gamma(getMaterialProperty<Real>("h_gamma")),
    // digonal terms of the first derivative of the interpolation function
    _dhalpha_dphialpha(getMaterialProperty<Real>("dhalpha_dphialpha")),
    _dhbeta_dphibeta(getMaterialProperty<Real>("dhbeta_dphibeta")),
    _dhgamma_dphigamma(getMaterialProperty<Real>("dhgamma_dphigamma")),
    
    // Non-digonal components of the interpolation matrix
    // for coupled variable phase_alpha
    _dhbeta_dphialpha(getMaterialProperty<Real>("dhbeta_dphialpha")),
    _dhgamma_dphialpha(getMaterialProperty<Real>("dhgamma_dphialpha")),
    
    // for coupled variable phase_beta
    _dhalpha_dphibeta(getMaterialProperty<Real>("dhalpha_dphibeta")),
    _dhgamma_dphibeta(getMaterialProperty<Real>("dhgamma_dphibeta")),
    
    // for coupled variable phase_gamma
    
    _dhalpha_dphigamma(getMaterialProperty<Real>("dhalpha_dphigamma")),
    _dhbeta_dphigamma(getMaterialProperty<Real>("dhbeta_dphigamma"))
    
{
}

Real
NMPhaseComposition::computeQpResidual()
{ 
  // Right-hand side of the equation
  
 const Real _weighted_sum = (_phase_comp_alpha[_qp] * _h_alpha[_qp] 
                     + _phase_comp_beta[_qp] * _h_beta[_qp] 
                     + _u[_qp] * _h_gamma[_qp]);

  // w_i*(sum(c_{theta}* h_theta} - c)= 0
  return _test[_i][_qp] * (_weighted_sum - _mole_fraction[_qp]);
}

Real
NMPhaseComposition::computeQpJacobian()
{
  return  _test[_i][_qp] * _h_gamma[_qp] * _phi[_j][_qp];
}

Real
NMPhaseComposition::computeQpOffDiagJacobian(unsigned int jvar)
{
  
  if (jvar == _phase_comp_beta_var)
  {  
    return _test[_i][_qp] * _h_beta[_qp] * _phi[_j][_qp];
  } 
  else if (jvar == _phase_comp_alpha_var)
  {
    return _test[_i][_qp] * _h_alpha[_qp] * _phi[_j][_qp];
  }  
  else if (jvar == _phase_alpha_var)
  {  
    return _test[_i][_qp] * ( _phase_comp_alpha[_qp] * _dhalpha_dphialpha[_qp] 
         + _phase_comp_beta[_qp] * _dhbeta_dphialpha[_qp] + _u[_qp] *_dhgamma_dphialpha[_qp]) * _phi[_j][_qp];
  }     
  else if (jvar == _phase_beta_var)
  {
    return _test[_i][_qp] *( _phase_comp_beta[_qp] * _dhbeta_dphibeta[_qp] 
        + _phase_comp_alpha[_qp] * _dhalpha_dphibeta[_qp] + _u[_qp] * _dhgamma_dphibeta[_qp]) * _phi[_j][_qp];
  }  
  else if (jvar == _phase_gamma_var)
  {
    return _test[_i][_qp] *(_u[_qp] * _dhgamma_dphigamma[_qp] 
        + _phase_comp_alpha[_qp] * _dhalpha_dphigamma[_qp] + _phase_comp_beta[_qp] *_dhbeta_dphigamma[_qp])* _phi[_j][_qp];
  }  
  else if (jvar == _mole_fraction_var)
  {
    return -_test[_i][_qp] * _phi[_j][_qp];
  }   
  else
  
    return 0.0;
}
