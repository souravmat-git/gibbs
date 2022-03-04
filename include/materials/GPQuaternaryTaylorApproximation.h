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
class GPQuaternaryTaylorApproximation;

//Included dependencies
#include "GPTernaryTaylorApproximation.h"

//template <>
//InputParameters validParams<GPQuaternaryTaylorApproximation>();

/**
 * Interpolation function is a Material class
 * that provides the free_energy of
 * of each phase alpha and beta
 */

class GPQuaternaryTaylorApproximation : public GPTernaryTaylorApproximation
{
public:
  GPQuaternaryTaylorApproximation(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:

   //Note: that we require more material property than KKS
   std::string _xD_name, _inv_D_tf_name;

   const Real & _xD_eqm;
   const Real & _D_tf_eqm;
   const Real & _D_diff_pot_eqm;

   MaterialProperty<Real> & _xD_val;
   MaterialProperty<Real> & _inv_D_tf_val;

  /**
   * Note: the coupled variable is mu
   */
  const VariableValue & _D_diff_pot;

};
