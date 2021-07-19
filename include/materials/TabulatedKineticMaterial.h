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

//#ifndef TABULATEDKINETICMATERIAL_H
//#define TABULATEDKINETICMATERIAL_H

#pragma once

// Forward Declarations
class TabulatedKineticMaterial;

//MOOSE includes
#include "Material.h"

template <>
InputParameters validParams<TabulatedKineticMaterial>();


class TabulatedKineticMaterial : public Material
{
public:
  TabulatedKineticMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override{};
    
  //Characteristic time:  to convert the data into non-dimensional terms
  const Real _tc; 
  
  //Characteristic length:  to convert the data into non-dimensional terms
  const Real _xc; 
  
  //Charctertistic energy J/mol
  const Real _Ec;
  
  //Molarvolume
  const Real _Vm;
   
};
//#endif // TABULATEDKINETICMATERIAL_H
