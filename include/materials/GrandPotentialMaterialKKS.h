//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef GRANDPOTENTIALMATERIALKKS_H
//#define GRANDPOTENTIALMATERIALKKS_H

#pragma once
// Forward Declarations
class GrandPotentialMaterialKKS;

//Included dependencies
#include "Material.h"

//template <>
//InputParameters validParams<GrandPotentialMaterialKKS>();

/**
 * Interpolation function is a Material class
 * that provides the free_energy of
 * of each phase alpha and beta
 */

class GrandPotentialMaterialKKS : public Material
{
public:
  GrandPotentialMaterialKKS(const InputParameters & parameters);
  
  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:
  /**
   * Define the member references that will hold the computed values
   * for the Real value properties in this class.
   * _h : holds the interpolation function
   * _dh : holds the first derivative of the interpolation function
   *_d2h : holds the second derivative of the interpolation function
   */

   //Note: that we require more material property than KKS

    MaterialProperty<Real> & _xB_alpha;
    MaterialProperty<Real> & _xB_beta;
    MaterialProperty<Real> & _inv_B_tf_alpha;
    MaterialProperty<Real> & _inv_B_tf_beta;
    MaterialProperty<Real> & _A_chem_pot_alpha;
    MaterialProperty<Real> & _A_chem_pot_beta;

  /**
   * Note: the coupled variable is mu
   */
  const VariableValue & _mu;

};
//#endif // GRANDPOTENTIALMATERIALKKS_H
