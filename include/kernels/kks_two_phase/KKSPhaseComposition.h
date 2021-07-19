//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef KKSPHASECOMPOSITION_H
#define KKSPHASECOMPOSITION_H

//Include dependencies
#include "Kernel.h"

// Forward Declarations
class KKSPhaseComposition;

template <>
InputParameters validParams<KKSPhaseComposition>();

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
 * \see KKSPhaseChemicalPotential
 * \see KKSHEtaPolyMaterial
 */
class KKSPhaseComposition : public Kernel
{
public:
  KKSPhaseComposition(const InputParameters & parameters);

protected:
  
  // The override command ensures that a virtual member function
  // from a base class is overridden
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:

  const VariableValue & _xB_alpha;
  unsigned int _xB_alpha_var;

  const VariableValue & _mole_fraction;
  unsigned int _mole_fraction_var;

  const VariableValue & _eta;
  unsigned int _eta_var;

  //_h is a reference variable that indicates the interpolation function
  const MaterialProperty<Real> & _h;
  // _dh is a reference variable that indicates the first derivative of the 
  // interpolation function
  const MaterialProperty<Real> &_dh;
  
};

#endif // KKSPHASECOMPOSITION_H
