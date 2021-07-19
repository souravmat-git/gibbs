//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This code is based on the code 
//*from TabulatedFluidProperties.C

//#ifndef TABULATEDPHASEMATERIAL_H
//#define TABULATEDPHASEMATERIAL_H

#pragma once

// Forward Declarations
class TabulatedPhaseMaterial;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<TabulatedPhaseMaterial>();


class TabulatedPhaseMaterial : public Material
{
public:
  TabulatedPhaseMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override{};
  
  //Molar volume to convert the thermodynamic data into J/m^{3} 
  const Real & _Vm;
    
  //Characteristic energy:  to convert the data into non-dimensional terms
  const Real & _Ec; 
   
};
//#endif // TABULATEDPHASEMATERIAL_H
