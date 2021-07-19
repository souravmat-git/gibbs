//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef KKSPHASEDIFFUSIONPOTENTIALC_H
#define KKSPHASEDIFFUSIONPOTENTIALC_H

#include "Kernel.h"


// Forward Declarations
class KKSPhaseDiffusionPotentialC;

template <>
InputParameters validParams<KKSPhaseDiffusionPotentialC>();

/**
 * Enforce the equality of the diffusion potentials in the two phases.
 * Eq. (21) in the original KKS paper.
 *
 * \f$ dF_beta/dc_beta - dF_alpha/dc_alpha =0\f$
 *
 */
class KKSPhaseDiffusionPotentialC : public Kernel
{
  public:
    KKSPhaseDiffusionPotentialC(const InputParameters & parameters);

  protected:
    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  private:
    /// coupled variable for xC_{beta}
    const  VariableValue & _xC_beta;
    unsigned int _xC_beta_var;
    
    // coupled variable for xB_beta
    const VariableValue & _xB_beta;
    unsigned int _xB_beta_var;
    
    //coupled variable for xB_alpha
    const VariableValue & _xB_alpha;
    unsigned int _xB_alpha_var;
    
    // coupled variable for xD_beta
    const VariableValue & _xD_beta;
    unsigned int _xD_beta_var;
    
    //coupled variable for xD_alpha
    const VariableValue & _xD_alpha;
    unsigned int _xD_alpha_var;
    
    /// material properties we need to access
    const MaterialProperty<Real> & _C_diff_pot_alpha;
    const MaterialProperty<Real> & _C_diff_pot_beta;
    
    const MaterialProperty<Real> & _C_therm_factor_alpha;
    const MaterialProperty<Real> & _C_therm_factor_beta;
    
    //To take into the cross terms in the equation
    //for ternary and quartenary alloys this kernel
    //needs to be modified
    const MaterialProperty<Real> & _BC_therm_factor_beta;
    const MaterialProperty<Real> & _BC_therm_factor_alpha;
    
    //Note that BD = DB 
    const MaterialProperty<Real> & _CD_therm_factor_beta;
    const MaterialProperty<Real> & _CD_therm_factor_alpha;
};
#endif // KKSPHASEDIFFUSIONPOTENTIALC_H
