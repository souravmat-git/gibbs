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

class StrainEnergyDrivingForce;

template<>
InputParameters validParams<StrainEnergyDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class StrainEnergyDrivingForce :  public Kernel
{
  public: 
    StrainEnergyDrivingForce(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
  //Displacement in the x-direction is the coupled variable
  //const VariableGradient & _grad_ux;
  unsigned int _ux_var;
  
  //Displacement in the y-direction is another coupled variable
  //const VariableGradient & _grad_uy;
  unsigned int _uy_var;
 
  //Material property required by the kernel 
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  //Elastic strain energy of the two phases
  const MaterialProperty<Real> & _fel_alpha;
  const MaterialProperty<Real> & _fel_beta;
  
  //Stress components sxx
  const MaterialProperty<Real> & _sxx_alpha;
  const MaterialProperty<Real> & _sxx_beta;
  
  //Stress components syy
  const MaterialProperty<Real> & _syy_alpha;
  const MaterialProperty<Real> & _syy_beta;
  
  //Stress components sxy
  const MaterialProperty<Real> & _sxy_alpha;
  const MaterialProperty<Real> & _sxy_beta; 
  
  //The phase-field mobility assumed to be constant
  const MaterialProperty<Real> & _L;
  
  //A non-dimensional factor 
  const MaterialProperty<Real> & _nd_factor;
};
