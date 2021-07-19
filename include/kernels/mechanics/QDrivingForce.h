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

class QDrivingForce;

template<>
InputParameters validParams<QDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class QDrivingForce :  public Kernel
{
  public: 
    QDrivingForce(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
  
  Real fdf() const;
  

  unsigned int _ux_var;
 
  //Material property required by the kernel
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  const MaterialProperty<Real> & _fel_alpha;
  const MaterialProperty<Real> & _fel_beta;
  
  const MaterialProperty<Real> & _sx_alpha;
  const MaterialProperty<Real> & _sx_beta;
  
  const MaterialProperty<Real> & _a;

  const MaterialProperty<Real> & _mat_const_alpha;
  const MaterialProperty<Real> & _mat_const_beta;
  
  const MaterialProperty<Real> & _L;
  
  //nd_factor
  const MaterialProperty<Real> & _nd_factor;
  
};
