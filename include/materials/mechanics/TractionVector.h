//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class TractionVector;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<TractionVector>();

class TractionVector : public Material
{
public:
  TractionVector(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    std::string _phase_name;
    
    //Independent variable on which this property depends
    const VariableGradient & _grad_eta;
    
   //Given stress components
    const MaterialProperty<Real> & _sxx;
    const MaterialProperty<Real> & _syy;
    const MaterialProperty<Real> & _sxy;
    
    std::string _tx_name, _ty_name;
    
    //Declare the traction vector
    MaterialProperty<Real> &  _tx;
    MaterialProperty<Real> &  _ty;  
};
