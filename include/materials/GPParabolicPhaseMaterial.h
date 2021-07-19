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
class GPParabolicPhaseMaterial;

//Included dependencies
#include "TabulatedPhaseMaterial.h"

template <>
InputParameters validParams<GPParabolicPhaseMaterial>();

/**
 * Interpolation function is a Material class
 * that provides the free_energy of 
 * of each phase alpha and beta 
 */
 
class GPParabolicPhaseMaterial : public TabulatedPhaseMaterial
{
public:
  GPParabolicPhaseMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
   
   //Note: that we require more material property than KKS

   std::string _xB_name, _inv_B_tf_name, _A_chem_pot_name;

   const Real & _A;
   const Real & _h;
   const Real & _k;
   
   MaterialProperty<Real> & _xB_val;
   MaterialProperty<Real> & _inv_B_tf_val;
   MaterialProperty<Real> & _A_chem_pot_val;
  /**
   * Note: the coupled variable is mu
   */
  const VariableValue & _B_diff_pot;

};
