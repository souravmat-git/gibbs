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

class KHSDrivingForce;

template<>
InputParameters validParams<KHSDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class KHSDrivingForce :  public Kernel
{
  public: 
    KHSDrivingForce(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
  //Displacement in the x-direction is the coupled variable
  const VariableGradient & _grad_ux;
  unsigned int _ux_var;
 
  //Material property required by the kernel 
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  const MaterialProperty<Real> & _elastic_energy_alpha;
  const MaterialProperty<Real> & _elastic_energy_beta;
  
  const MaterialProperty<Real> & _sx_alpha;
  const MaterialProperty<Real> & _sx_beta;
  
  const MaterialProperty<Real> & _L;
  
  //nd_factor
  const MaterialProperty<Real> & _nd_factor;
  
  
};
