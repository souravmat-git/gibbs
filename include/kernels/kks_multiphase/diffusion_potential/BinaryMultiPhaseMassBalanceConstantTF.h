//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef BINARYMULTIPHASEMASSBALANCE_H
//#define BINARYMULTIPHASEMASSBALANCE_H

#pragma once

#include "BinaryMultiPhaseMassBalance.h"

class BinaryMultiPhaseMassBalanceConstantTF;

template <>
InputParameters validParams<BinaryMultiPhaseMassBalanceConstantTF>();

/**
 * Kernel to implement the continuity
 * equation for mass conservation
 * variable on which this kernel operates: X
 */
class BinaryMultiPhaseMassBalanceConstantTF: public BinaryMultiPhaseMassBalance
{
public:
  BinaryMultiPhaseMassBalanceConstantTF(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;


};
