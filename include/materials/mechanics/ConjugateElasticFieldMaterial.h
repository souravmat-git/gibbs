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
class ConjugateElasticFieldMaterial;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<ConjugateElasticFieldMaterial>();

class ConjugateElasticFieldMaterial : public Material
{
public:
  ConjugateElasticFieldMaterial(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the strain, stress and the elastic energy
  //of a phase with or without a prescribed eigenstrain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Independent variable on which this property depends
    const VariableValue & _stress_x;

    //Input material properties
    //for this material class are material constant and transformation strain
    const MaterialProperty<Real> & _eT;
    const MaterialProperty<Real> & _mat_const;

    //Complementary stress energy, its first derivative is strain and compliance
    MaterialProperty<Real> & _ex_val;
    MaterialProperty<Real> & _compliance_val;
    MaterialProperty<Real> & _comp_energy_val;

};
