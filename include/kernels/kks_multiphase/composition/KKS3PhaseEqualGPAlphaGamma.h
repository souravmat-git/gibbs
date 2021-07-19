//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef KKS3PHASEEQUALGPALPHAGAMMA_H
#define KKS3PHASEEQUALGPALPHAGAMMA_H

#include "Kernel.h"

class KKS3PhaseEqualGPAlphaGamma;

template<>
InputParameters validParams< KKS3PhaseEqualGPAlphaGamma>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh_{alpha}/dphi_{gamma}} * (omega_alpha - omega_gamma)
  *omega_beta - Grandpotential of the gamma phase
  *omega_alpha - GrandPotential of the alpha phase
  **/

class  KKS3PhaseEqualGPAlphaGamma:  public Kernel
{
  public: 
     KKS3PhaseEqualGPAlphaGamma(const InputParameters & parameters);
  
  protected:
  
    // The override command ensures that a virtual member function
    // from a base class is overridden
    
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

//private: Member function defined within private cannot be used 
// for derived class

    const VariableValue & _phase_comp_alpha;
    unsigned int _phase_comp_alpha_var;
    
    const VariableValue & _phase_comp_gamma;
    unsigned int _phase_comp_gamma_var;
    
    const VariableValue & _phase_beta;
    unsigned int _phase_beta_var;
    
    const VariableValue & _phase_alpha;
    unsigned int _phase_alpha_var;
        
    const MaterialProperty<Real> & _free_energy_alpha;
    const MaterialProperty<Real> & _free_energy_gamma;
    const MaterialProperty<Real> & _df_alpha;
    const MaterialProperty<Real> & _df_gamma;
    const MaterialProperty<Real> & _d2f_alpha;
    const MaterialProperty<Real> & _d2f_gamma;

    
    //In addition this kernel requires the first and second derivatives
    // of the interpolation function
    const MaterialProperty<Real> & _dhalpha_dphigamma;
    const MaterialProperty<Real> & _d2halpha_dphigamma2;
    const MaterialProperty<Real> & _d2halpha_dphigamma_dphibeta;
    const MaterialProperty<Real> & _d2halpha_dphigamma_dphialpha;
    
    //The mobility here is assumed to be a constant
    const MaterialProperty<Real> & _L;
    
};
#endif // KKS3PHASEEQUALGPALPHAGAMMA_H
