//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class RelativeDisplacementPerUnitLength;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<RelativeDisplacementPerUnitLength>();

class RelativeDisplacementPerUnitLength : public Material
{
public:
  RelativeDisplacementPerUnitLength(const InputParameters & parameters);

protected:
  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    std::string _phase_name;

    //Independent variable on which this property depends
    const VariableGradient & _grad_eta;
    
    //Obtain the displacement gradient matrix from the displacement vector
    const VariableGradient & _grad_ux;
    const VariableGradient & _grad_uy;
    
    std::string _dlx_name, _dly_name;
      
    //Declare the deformation vector
    MaterialProperty<Real> &  _dlx;
    MaterialProperty<Real> &  _dly;  
};
