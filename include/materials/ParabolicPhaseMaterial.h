//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef PARABOLICPHASEMATERIAL_H
//#define PARABOLICPHASEMATERIAL_H

#pragma once
// Forward Declarations
class ParabolicPhaseMaterial;

//Included dependencies
#include "Material.h"

template <>
InputParameters validParams<ParabolicPhaseMaterial>();

/* A material class to create a vertex form of a parabola
 * y = 0.5*A*(x-h)^{2} + k
 * A = curvature of the parabola
 * h = horizontal shift of the parabola
 * k = vertical shift of the parabola
 * to supply the free energy, diffusion potential & therm_factor
 * for a binary alloy A-B, where B is the independent comp
 */

 
class ParabolicPhaseMaterial: public Material
{
public:
  ParabolicPhaseMaterial(const InputParameters & parameters);

protected:

  virtual void computeQpProperties() override;
  
private:

  //String variable to hold th free energy, diffusion potential
  //and the thermodynamic facor name obtained from the input file
  std::string _f_phase_name, _B_diff_pot_phase_name, _B_therm_factor_phase_name;
    
  MaterialProperty<Real> & _f_phase_val;
  MaterialProperty<Real> & _B_diff_pot_phase_val;
  MaterialProperty<Real> & _B_therm_factor_phase_val;
  
  //Independent variable on which this property depends
  const VariableValue & _mol_fraction_B;
  
  //This class needs three parameters
  //The curvature of the parabola
  const Real _A;
  //The horizontal shift of the parabola
  const Real _h;
  //The vertical shift of the parabola
  const Real _k;
  
};
//#endif // PARABOLICPHASEMATERIAL_H
