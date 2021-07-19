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

class SpecifyTangentialDisplacementBC;

template <>
InputParameters validParams<SpecifyTangentialDisplacementBC>();

/**
 * Implements the radial displacement which varies as a function of 
 */
class SpecifyTangentialDisplacementBC : public NodalBC
{
public:
  SpecifyTangentialDisplacementBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian (unsigned int jvar) override;

private:

  /// x -displacement
  const VariableValue & _ux;  
  
  //To couple the y-displacement
  unsigned int _ux_var;
  
  //Real tangential displacement
  const Real & _uT;
  
  //Center point
  const Point _center_point;
  
};

