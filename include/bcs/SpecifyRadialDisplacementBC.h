//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "NodalBC.h"

class SpecifyRadialDisplacementBC;

template <>
InputParameters validParams<SpecifyRadialDisplacementBC>();

/**
 * Implements the radial displacement which varies as a function of 
 */
class SpecifyRadialDisplacementBC : public NodalBC
{
public:
  SpecifyRadialDisplacementBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian (unsigned int jvar) override;

private:

  /// y -displacement
  const VariableValue & _uy;  
  
  //To couple the y-displacement
  unsigned int _uy_var;
  
  //Real radial displacement
  const Real & _uR;
  
  //Center point
  const Point _center_point;
  
};

