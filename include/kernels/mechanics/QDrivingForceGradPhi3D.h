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

class QDrivingForceGradPhi3D;

template<>
InputParameters validParams<QDrivingForceGradPhi3D>();

/** Implements the gradient driving force 
  **/

class QDrivingForceGradPhi3D :  public Kernel
{
  public: 
    QDrivingForceGradPhi3D(const InputParameters & parameters);
  
  protected:
  
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
  
  //df_dgrad_phi and d2f_dgradphi
  const MaterialProperty<RealVectorValue> & _df_dgradphi;
  const MaterialProperty<RealVectorValue> & _d2f_dgradphi_dphi;
  
  //Mobility
  const MaterialProperty<Real> & _L;
  
  //nd_factor
  const MaterialProperty<Real> & _nd_factor;
};
