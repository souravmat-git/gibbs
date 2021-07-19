//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ALPHAGAMMAEQUALDIFFUSIONPOTENTIAL_H
#define ALPHAGAMMAEQUALDIFFUSIONPOTENTIAL_H

#include "Kernel.h"

// Forward Declarations
class AlphaGammaEqualDiffusionPotential;

template <>
InputParameters validParams<AlphaGammaEqualDiffusionPotential>();

/**
 * Enforce the equality of the diffusion potentials in the two phases.
 * Eq. (21) in the original KKS paper.
 *
 * \f$ dF_beta/dc_beta - dF_gamma/dc_gamma =0\f$
 *
 */
class AlphaGammaEqualDiffusionPotential : public Kernel
{
  public:
    AlphaGammaEqualDiffusionPotential(const InputParameters & parameters);

  protected:
    virtual Real computeQpResidual();
    virtual Real computeQpJacobian();
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  private:
    /// coupled variable for c_{gamma}
    const  VariableValue & _phase_comp_gamma;
    unsigned int _phase_comp_gamma_var;

    /// material properties we need to access
    const MaterialProperty<Real> & _df_gamma;
    const MaterialProperty<Real> & _df_alpha;
    const MaterialProperty<Real> & _d2f_gamma;
    const MaterialProperty<Real> & _d2f_alpha;
};
#endif //ALPHAGAMMAEQUALDIFFUSIONPOTENTIAL_H
