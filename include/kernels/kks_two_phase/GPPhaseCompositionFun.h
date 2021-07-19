//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//Include dependencies
#include "Kernel.h"

// Forward Declarations
class GPPhaseCompositionFun;
class Function;

template <>
InputParameters validParams<GPPhaseCompositionFun>();

/**
 * Enforce sum of phase concentrations to be the real concentration.
 *
 * \f$ c=h(\phi)c_{beta}+\left(1-h(\phi)\right)c_{alpha}\f$
 *
 * The non-linear variable for this Kernel is the concentration \f$ c_{beta} \f$, while
 * \f$ c_alpha \f$ and \f$ c \f$ are supplied as coupled variables.
 * (compare this to KKSPhaseChemicalPotential, where the non-linear variable is
 * the other phase concentration \f$ c_a \f$!)
 *
 */
class GPPhaseCompositionFun : public Kernel
{
public:
  GPPhaseCompositionFun(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  
  //Two material properties i.e the phase-mole fractions
  const MaterialProperty<Real> & _xB_alpha;
  const MaterialProperty<Real> & _xB_beta;
  
  const MaterialProperty<Real> & _inv_B_tf_alpha;
  const MaterialProperty<Real> & _inv_B_tf_beta;
  
  const Function & _xB;
  const Function & _eta;

  const MaterialProperty<Real> & _h;

};
