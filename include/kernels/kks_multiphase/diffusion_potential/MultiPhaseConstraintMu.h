//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef MULTIPHASECONSTRAINTMU_H
//#define MULTIPHASECONSTRAINTMU_H

#pragma once

#include "Kernel.h"

// Forward declaration
class MultiPhaseConstraintMu;

//class template specialization:
template <>
InputParameters validParams<MultiPhaseConstraintMu>();

/**
 * This class enforces the equation for mole fraction
 * c = c_{alpha}* h_alpha + c_{beta}* h_beta + c_{gamma} * h_gamma
 * The class is a derived class from base class kernel
 **/
class MultiPhaseConstraintMu : public Kernel
{
public:
   MultiPhaseConstraintMu(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;
       
     
  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;
  
  const VariableValue & _phase_gamma;
  unsigned int _phase_gamma_var;
  
  const VariableValue & _phase_delta;
  unsigned int _phase_delta_var;
  
  const VariableValue & _phase_epsilon;
  unsigned int _phase_epsilon_var;
          
  const MaterialProperty<Real> & _h_alpha;
  const MaterialProperty<Real> & _h_beta;
  const MaterialProperty<Real> & _h_gamma;
  const MaterialProperty<Real> & _h_delta;
  const MaterialProperty<Real> & _h_epsilon;    
                
  const MaterialProperty<Real> & _dhbeta_dphialpha;
  const MaterialProperty<Real> & _dhgamma_dphialpha;
  const MaterialProperty<Real> & _dhdelta_dphialpha;
  const MaterialProperty<Real> & _dhepsilon_dphialpha; 
        
  const MaterialProperty<Real> & _dhalpha_dphibeta;
  const MaterialProperty<Real> & _dhgamma_dphibeta;
  const MaterialProperty<Real> & _dhdelta_dphibeta;
  const MaterialProperty<Real> & _dhepsilon_dphibeta;
        
  const MaterialProperty<Real> & _dhalpha_dphigamma;
  const MaterialProperty<Real> & _dhbeta_dphigamma;
  const MaterialProperty<Real> & _dhdelta_dphigamma;
  const MaterialProperty<Real> & _dhepsilon_dphigamma;
  
  const MaterialProperty<Real> & _dhalpha_dphidelta;
  const MaterialProperty<Real> & _dhbeta_dphidelta;
  const MaterialProperty<Real> & _dhgamma_dphidelta;
  const MaterialProperty<Real> & _dhepsilon_dphidelta;
  
  const MaterialProperty<Real> & _dhalpha_dphiepsilon;
  const MaterialProperty<Real> & _dhbeta_dphiepsilon;
  const MaterialProperty<Real> & _dhgamma_dphiepsilon;
  const MaterialProperty<Real> & _dhdelta_dphiepsilon; 
  
  const VariableValue & _xB;
  unsigned int _xB_var;
  
  //phase comp material property for each phase
  const MaterialProperty<Real> & _xB_alpha;
  const MaterialProperty<Real> & _xB_beta;
  const MaterialProperty<Real> & _xB_gamma;
  const MaterialProperty<Real> & _xB_delta;
  const MaterialProperty<Real> & _xB_epsilon;
   
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> & _inv_B_tf_beta;
  const MaterialProperty<Real> & _inv_B_tf_gamma;
  const MaterialProperty<Real> & _inv_B_tf_delta;
  const MaterialProperty<Real> & _inv_B_tf_epsilon;
    
};
//#endif // MULTIPHASECONSTRAINTMU_H
