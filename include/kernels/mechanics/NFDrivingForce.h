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

class NFDrivingForce;

template<>
InputParameters validParams<NFDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class NFDrivingForce :  public Kernel
{
  public: 
     NFDrivingForce(const InputParameters & parameters);
  
  protected:
  
  Real QuantMechDrivingForce() const;
  Real dR_dex_alpha() const;
  Real dR_dex_beta() const;
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
  //Strain the x-direction is the coupled variable
  const VariableValue & _ex_alpha;
  unsigned int _ex_alpha_var;
  
  const VariableValue & _ex_beta;
  unsigned int _ex_beta_var;
 
  //Material property required by the kernel 
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  const MaterialProperty<Real> & _fel_alpha;
  const MaterialProperty<Real> & _fel_beta;
  const MaterialProperty<Real> & _sx_alpha;
  const MaterialProperty<Real> & _sx_beta;
  const MaterialProperty<Real> & _mat_const_alpha;
  const MaterialProperty<Real> & _mat_const_beta;
  
  const MaterialProperty<Real> & _L;
  
  const MaterialProperty<Real> & _nd_factor;
};
