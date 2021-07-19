//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "AuxKernel.h"


/**
 * Calculates the radial displacement from the cartesian 
 * components. However, it is limted to polar coordinates.
 */

class CylindricalRankOneAux : public AuxKernel
{
public:
  static InputParameters validParams();

  CylindricalRankOneAux(const InputParameters & parameters);
  virtual ~CylindricalRankOneAux() {};

protected:
  virtual Real computeValue() override;
  
  const VariableValue & _disp_x;
  const VariableValue & _disp_y;
  
  //variable to store the specific component
  const unsigned int _component;
  
  //Center point location
  const Point _center_point;
};
