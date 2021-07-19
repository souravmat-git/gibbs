//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class DeformationVector;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<DeformationVector>();

class DeformationVector : public Material
{
public:
  DeformationVector(const InputParameters & parameters);

protected:
  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    std::string _phase_name;


    //Independent variable on which this property depends
    const VariableGradient & _grad_eta;
    
   //Given stress components
    const MaterialProperty<Real> & _exx;
    const MaterialProperty<Real> & _eyy;
    const MaterialProperty<Real> & _exy;
    
    std::string _ax_name, _ay_name;
    
    //Declare the deformation vector
    MaterialProperty<Real> &  _ax;
    MaterialProperty<Real> &  _ay;  
};
