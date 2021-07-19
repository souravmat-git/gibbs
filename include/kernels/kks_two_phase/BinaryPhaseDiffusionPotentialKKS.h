//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BINARYPHASEDIFFUSIONPOTENTIALKKS_H
#define BINARYPHASEDIFFUSIONPOTENTIALKKS_H

#include "Kernel.h"


// Forward Declarations
class BinaryPhaseDiffusionPotentialKKS;

template <>
InputParameters validParams<BinaryPhaseDiffusionPotentialKKS>();

/**
 * Enforce the equality of diffusion potentials for component B in the two phases.
 * Eq. (21) in the original KKS paper.
 * \tilde{\mu}_{B}^{beta} - \tilde{\mu}_{B}^{alpha} = 0
 * can be used for binary alloys, ternary and quaternary alloys 
 */
class BinaryPhaseDiffusionPotentialKKS : public Kernel
{
  public:
    BinaryPhaseDiffusionPotentialKKS(const InputParameters & parameters);

  protected:
  
    virtual Real computeQpResidual() override;
    virtual Real computeQpJacobian() override;
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    /// coupled variable for xB_{beta}
    const  VariableValue & _xB_beta;
    unsigned int _xB_beta_var;
  
    /// Diffusion potentials of component B in each phase
    const MaterialProperty<Real> & _B_diff_pot_alpha;
    const MaterialProperty<Real> & _B_diff_pot_beta;
    
    //Therm. factor of comp B in each phase
    const MaterialProperty<Real> & _B_therm_factor_alpha;
    const MaterialProperty<Real> & _B_therm_factor_beta;
    
};
#endif // BINARYPHASEDIFFUSIONPOTENTIALKKS_H
