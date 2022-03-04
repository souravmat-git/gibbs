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
class InterpolationEtaFunction;
class Function;

//Included dependencies
#include "Material.h"

//template <>
//InputParameters validParams<InterpolationEtaFunction>();

/**
 * Interpolation function is a Material class
 * that provides the interpolation function required by the kernel
 * the form of the interpolation function is
 * h(eta) = eta^(3)(6 * eta^(2) - 15*eta + 10.0)
 */

class InterpolationEtaFunction : public Material
{
public:
  InterpolationEtaFunction(const InputParameters & parameters);

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
   MaterialProperty<Real> & _h;
   MaterialProperty<Real> & _dh;
   MaterialProperty<Real> & _d2h;

   //MaterialProperty<Real> & _g;
  /**
   * This is the member reference that will hold the coupled variable
   */
  const Function & _eta;
};
