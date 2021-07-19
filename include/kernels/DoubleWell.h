//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DOUBLEWELL_H
#define DOUBLEWELL_H

#include "Kernel.h"

class DoubleWell;

template <>
InputParameters validParams<DoubleWell>();

/**
 * Kernel to implement double well
 */
class DoubleWell : public Kernel
{
public:
  DoubleWell(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Mobility & Barrier height
  const MaterialProperty<Real> & _BH;
  const MaterialProperty<Real> & _L;
};

#endif // DOUBLEWELL_H
