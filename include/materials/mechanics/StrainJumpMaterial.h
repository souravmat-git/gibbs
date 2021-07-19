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
class StrainJumpMaterial;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<StrainJumpMaterial>();

class StrainJumpMaterial : public Material
{
public:
  StrainJumpMaterial(const InputParameters & parameters);
  
  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    //Independent variable on which this property depends
    const VariableGradient & _grad_ux;
    
    //Input material properties
    const MaterialProperty<Real> & _mat_const_alpha;
    const MaterialProperty<Real> & _mat_const_beta;
    const MaterialProperty<Real> & _eT_alpha;
    const MaterialProperty<Real> & _eT_beta;
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _dh;
    
    MaterialProperty<Real> & _a_val;
    MaterialProperty<Real> & _da_de_val;
    MaterialProperty<Real> & _da_dh_val;  
};
