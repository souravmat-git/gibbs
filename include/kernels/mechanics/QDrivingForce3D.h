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

class QDrivingForce3D;

template<>
InputParameters validParams<QDrivingForce3D>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class QDrivingForce3D :  public Kernel
{
  public: 
    QDrivingForce3D(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
  
  Real driving_force() const;
  
  //unsigned int _ux_var;
 
  //Material property required by the kernel
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  const MaterialProperty<Real> & _alpha_strain_energy;
  const MaterialProperty<Real> & _beta_strain_energy;
  
  const MaterialProperty<RankTwoTensor> & _alpha_stress;
  const MaterialProperty<RankTwoTensor> & _beta_stress;
  
  const MaterialProperty<RankTwoTensor> & _strain_jump;

  const MaterialProperty<RankFourTensor> & _alpha_stiffness;
  const MaterialProperty<RankFourTensor> & _beta_stiffness;
  
  const MaterialProperty<Real> & _L;
  
  //nd_factor
  const MaterialProperty<Real> & _nd_factor;
  
};
