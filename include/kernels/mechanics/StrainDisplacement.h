//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class StrainDisplacement;

//MOOSE includes
#include "Kernel.h"

template <>
InputParameters validParams<StrainDisplacement>();

class StrainDisplacement : public Kernel
{
public:
  StrainDisplacement(const InputParameters & parameters);
  
  const VariableValue & _ex;
  unsigned int _ex_var;
 

protected:

   virtual Real computeQpResidual() override;
   virtual Real computeQpJacobian() override;
   virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
   
};
