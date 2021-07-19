//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef KKSTWOPHASECOMPOSITION_H
#define KKSTWOPHASECOMPOSITION_H

//Include dependencies
#include "Kernel.h"

// Forward Declarations
class KKSTwoPhaseComposition;

template <>
InputParameters validParams<KKSTwoPhaseComposition>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 *
 * \f$ c=h_{beta} * c_{beta}+ h_alpha * c_{alpha}\f$
 *
 * The non-linear variable for this Kernel is the c_beta
 *

 */
class KKSTwoPhaseComposition : public Kernel
{
public:
  KKSTwoPhaseComposition(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const VariableValue & _phase_comp_alpha;
  unsigned int _phase_comp_alpha_var;

  const VariableValue & _mole_fraction;
  unsigned int _mole_fraction_var;

  const VariableValue & _phase_alpha;
  unsigned int _phase_alpha_var;
  
  const VariableValue & _phase_beta;
  unsigned int _phase_beta_var;
};

#endif // KKSTWOPHASECOMPOSITION_H
