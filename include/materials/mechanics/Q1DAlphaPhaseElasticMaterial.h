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
class Q1DAlphaPhaseElasticMaterial;
class Function;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<Q1DAlphaPhaseElasticMaterial>();

class Q1DAlphaPhaseElasticMaterial : public Material
{
public:
  Q1DAlphaPhaseElasticMaterial(const InputParameters & parameters);
  
  //This class calculates the strain, stress and the elastic energy
  //of a phase with or without a prescribed eigenstrain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //Independent variable on which the material property depends
    const VariableGradient & _grad_ux;

    //Input material properties
    //for this material class are mat_Const, transformation strain and strain jump
    const MaterialProperty<Real> & _a;
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _eT;
    const MaterialProperty<Real> & _mat_const;
    
    MaterialProperty<Real> & _ex_val;
    MaterialProperty<Real> & _sx_val;
    MaterialProperty<Real> & _elastic_energy_val;    
};
