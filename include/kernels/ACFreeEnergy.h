//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ACFREEENERGY_H
#define ACFREEENERGY_H

#include "Kernel.h"

class ACFreeEnergy;

template <>
InputParameters validParams<ACFreeEnergy>();

/**
 * Compute the Allen-Cahn interface term with constant Mobility and Interfacial parameter
 */
class ACFreeEnergy : public Kernel
{
public:
  ACFreeEnergy(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  /// Mobility
  const MaterialProperty<Real> & _L;
};

#endif // ACFREEENERGY_H
