//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class ObtainUnitNormalBase;

//MOOSE includes
#include "Kernel.h"

template <>
InputParameters validParams<ObtainUnitNormalBase>();

class ObtainUnitNormalBase : public Kernel
{
public:
  ObtainUnitNormalBase(const InputParameters & parameters);
    
protected:

   virtual Real computeQpResidual() override;
   virtual Real computeQpJacobian() override;
   virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
   
   Real nx()  const;
   Real ny()  const;
   Real nz()  const; 
  
   const Real _tol;
   const Real _val;
    
   //Obtain the gradient of phase-field variable
   const VariableGradient & _grad_eta;
   
};
