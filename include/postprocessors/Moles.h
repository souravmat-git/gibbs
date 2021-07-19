//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "ElementIntegralVariablePostprocessor.h"

// Forward Declarations
class Moles;

template <>
InputParameters validParams<Moles>();

/**
 * This postprocessor computes the total number of moles of a component
 * This class is derived from ElementVariablePostprocessor 
 * and it overrides the computeqpIntegral in the base class
 */
class Moles : public ElementIntegralVariablePostprocessor
{
public:
  Moles(const InputParameters & parameters);

protected:

  virtual Real computeQpIntegral() override;
  
  Real _volume; 
}; 

