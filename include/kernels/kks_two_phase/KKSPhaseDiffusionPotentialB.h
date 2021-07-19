//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef KKSPHASEDIFFUSIONPOTENTIALB_H
#define KKSPHASEDIFFUSIONPOTENTIALB_H

#include "BinaryPhaseDiffusionPotentialKKS.h"

// Forward Declarations
class KKSPhaseDiffusionPotentialB;

template <>
InputParameters validParams<KKSPhaseDiffusionPotentialB>();

/**
 * Enforce the equality of diffusion potentials for component B in the two phases.
 * Eq. (21) in the original KKS paper.
 * \tilde{\mu}_{B}^{beta} - \tilde{\mu}_{B}^{alpha} = 0
 * can be used for binary alloys, ternary and quaternary alloys 
 */
class KKSPhaseDiffusionPotentialB : public BinaryPhaseDiffusionPotentialKKS
{
  public:
    KKSPhaseDiffusionPotentialB(const InputParameters & parameters);

  protected:
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  private:
    
    //coupledVariable Component C in beta phase
    const VariableValue & _xC_beta;
    unsigned int _xC_beta_var;
    
    //coupledVariable Component C in alpha phase
    const VariableValue & _xC_alpha;
    unsigned int _xC_alpha_var;
    
    //coupledVariable Component D in beta phase
    const VariableValue & _xD_beta;
    unsigned int _xD_beta_var;
    
    //coupledVariable Component D in alpha phase
    const VariableValue & _xD_alpha;
    unsigned int _xD_alpha_var;

    //Thermfactor of component B,C mu(B).x(C)
    const MaterialProperty<Real> & _BC_therm_factor_beta;
    const MaterialProperty<Real> & _BC_therm_factor_alpha;
    
   //Thermfactor of component B,D mu(B).x(D)
    const MaterialProperty<Real> & _BD_therm_factor_beta;
    const MaterialProperty<Real> & _BD_therm_factor_alpha;
    
};
#endif // KKSPHASEDIFFUSIONPOTENTIALB_H
