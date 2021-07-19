//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* This code modifies the KKSPhaseConcentration in MOOSE
//* by including the hand coded interpolation function
//* instead of using the Material property- S.Chatterjee

#include "KKSTwoPhaseComposition.h"

registerMooseObject("gibbsApp", KKSTwoPhaseComposition);

template <>
InputParameters
validParams<KKSTwoPhaseComposition>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Eqn: (h_alpha)*c_alpha + h_beta*c_beta - c = 0."
                             "non-linear variable of this kernel is c_beta.");
  params.addRequiredCoupledVar("phase_comp_alpha", "Phase concentration in alpha phase");
  params.addRequiredCoupledVar("mole_fraction", "Real concentration");
  params.addRequiredCoupledVar("phase_alpha", "To distinguish phases");
  params.addRequiredCoupledVar("phase_beta","For beta phase");
  return params;
}

KKSTwoPhaseComposition::KKSTwoPhaseComposition(const InputParameters & parameters)
  :  Kernel(parameters),
    _phase_comp_alpha(coupledValue("phase_comp_alpha")),
    _phase_comp_alpha_var(coupled("phase_comp_alpha")),
    _mole_fraction(coupledValue("mole_fraction")),
    _mole_fraction_var(coupled("mole_fraction")),
    _phase_alpha(coupledValue("phase_alpha")),
    _phase_alpha_var(coupled("phase_alpha")),
    _phase_beta(coupledValue("phase_beta")),
    _phase_beta_var(coupled("phase_beta"))
{
}

Real
KKSTwoPhaseComposition::computeQpResidual()
{
  // R: w*((h_beta*c_beta + h_alpha* c_alpha -c)
  
  const Real _sum = (_phase_alpha[_qp] * _phase_alpha[_qp] 
                    + _phase_beta[_qp] * _phase_beta[_qp]);
                    
  //Interpolation function                    
                     
  const Real _h_alpha = std::pow(_phase_alpha[_qp],2.0)/_sum;  
  const Real _h_beta = std::pow(_phase_beta[_qp],2.0)/_sum;
  
  
  return _test[_i][_qp] * (_h_beta*_u[_qp] + _h_alpha * _phase_comp_alpha[_qp] 
         - _mole_fraction[_qp]);
}

Real
KKSTwoPhaseComposition::computeQpJacobian()
{
  const Real _sum = (_phase_alpha[_qp] * _phase_alpha[_qp] 
                    + _phase_beta[_qp] * _phase_beta[_qp]);
                    
  //Interpolation function for phase beta
  const Real _h_beta = std::pow(_phase_beta[_qp],2.0)/_sum;
   
  return _test[_i][_qp] * _h_beta * _phi[_j][_qp];
}

Real
KKSTwoPhaseComposition::computeQpOffDiagJacobian(unsigned int jvar)
{

  const Real _sum = (_phase_alpha[_qp] * _phase_alpha[_qp] 
                    + _phase_beta[_qp] * _phase_beta[_qp]);
                    
  //Interpolation function                    
                     
  const Real _h_alpha = std::pow(_phase_alpha[_qp],2.0)/_sum;  
  const Real _h_beta = std::pow(_phase_beta[_qp],2.0)/_sum;


  if (jvar == _phase_comp_alpha_var)
  {
    return _test[_i][_qp] * _h_alpha * _phi[_j][_qp];
  }
  else if (jvar == _mole_fraction_var)
  {
    return -_test[_i][_qp] * _phi[_j][_qp];
  }
  else if (jvar == _phase_alpha_var)
  {
    return (2 * _phase_alpha[_qp] /_sum) * (_h_beta * _u[_qp] + (1-_h_alpha) * _phase_comp_alpha[_qp]);
  }  
  else if (jvar == _phase_beta_var)
  {
    return (2 * _phase_beta[_qp]/_sum) * ((1-_h_beta) * _u[_qp] + _h_alpha * _phase_comp_alpha[_qp]);
  }
  else
    return 0.0;
}
