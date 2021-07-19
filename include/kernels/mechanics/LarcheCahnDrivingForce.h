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

class LarcheCahnDrivingForce;

template<>
InputParameters validParams<LarcheCahnDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the LC model
  *h^{\prime}(fc_beta - fc_alpha)= 0
  **/

class LarcheCahnDrivingForce :  public Kernel
{
  public: 
    LarcheCahnDrivingForce(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
  //Stress in the x-direction is a copled variable
  const VariableValue & _stress_x;
  unsigned int _stress_x_var;
 
  //Material property required by the kernel 
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  const MaterialProperty<Real> & _comp_energy_alpha;
  const MaterialProperty<Real> & _comp_energy_beta;
  
  const MaterialProperty<Real> & _ex_alpha;
  const MaterialProperty<Real> & _ex_beta;
  
  const MaterialProperty<Real> & _L;
  const MaterialProperty<Real> & _nd_factor;
};
