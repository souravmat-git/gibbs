//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class GPTernaryTaylorApproximation;

//Included dependencies
#include "GPTaylorApproximation.h"

//template <>
//InputParameters validParams<GPTernaryTaylorApproximation>();

/**
 * Interpolation function is a Material class
 * that provides the free_energy of
 * of each phase alpha and beta
 */

class GPTernaryTaylorApproximation : public GPTaylorApproximation
{
public:
  GPTernaryTaylorApproximation(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

   //Note: that we require more material property than KKS
   std::string _xC_name, _inv_C_tf_name;

   const Real & _xC_eqm;
   const Real & _C_tf_eqm;
   const Real & _C_diff_pot_eqm;

   MaterialProperty<Real> & _xC_val;
   MaterialProperty<Real> & _inv_C_tf_val;
  /**
   * Note: the coupled variable is mu
   */
  const VariableValue & _C_diff_pot;

};
