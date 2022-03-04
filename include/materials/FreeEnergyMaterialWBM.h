//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef FREEENERGYMATERIALWBM_H
//#define FREEENERGYMATERIALWBM_H

#pragma once

// Forward Declarations
class FreeEnergyMaterialWBM;

//Included dependencies
#include "Material.h"

//template <>
//InputParameters validParams<FreeEnergyMaterialWBM>();

/**
 * Interpolation function is a Material class
 * that provides the free_energy of
 * of each phase alpha and beta
 */

class FreeEnergyMaterialWBM : public Material
{
public:
  FreeEnergyMaterialWBM(const InputParameters & parameters);

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
    MaterialProperty<Real> & _free_energy_alpha;
    MaterialProperty<Real> & _dfalpha_dc;
    MaterialProperty<Real> & _d2falpha_dc2;

    MaterialProperty<Real> & _free_energy_beta;
    MaterialProperty<Real> & _dfbeta_dc;
    MaterialProperty<Real> & _d2fbeta_dc2;

  /**
   * Note that WBM has the same composition
   */
  const VariableValue & _comp;
};

//#endif // FREEENERGYMATERIALWBM_H
