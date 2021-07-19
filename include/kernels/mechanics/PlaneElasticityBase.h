//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Kernel.h"

class PlaneElasticityBase;

template <>
InputParameters validParams<PlaneElasticityBase>();

/**
 * Kernel to implement the momentum balance
 * in x-direction for a single phase material
 * The variable that this kernel acts on is the displacement field in the x-direction.
 * Further, this variable is coupled to the displacement in the y-direction
 */
class PlaneElasticityBase : public Kernel
{
public:
  PlaneElasticityBase(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  //C11, C12, C22, C66 
  const MaterialProperty<Real> & _C11;
  const MaterialProperty<Real> & _C12;
  const MaterialProperty<Real> & _C22;
  const MaterialProperty<Real> & _C66;
   
};
