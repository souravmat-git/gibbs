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

class PlaneElasticityTwoPhaseBase;

template <>
InputParameters validParams<PlaneElasticityTwoPhaseBase>();

/**
 * Kernel to precompute
 * the a) overall C11, C12, C22, C66
 * the b) and the eigenstress
 */
class PlaneElasticityTwoPhaseBase : public Kernel
{
public:
  PlaneElasticityTwoPhaseBase(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  const VariableValue & _eta;
  unsigned int _eta_var;
  
  std::string _h_name;

  //Obtain the interpolation function and its derivative
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  //Overall elastic constants
  Real C11_interp() const;
  Real C12_interp() const;
  Real C22_interp() const;
  Real C66_interp() const;
   
  
  //plane elasticity stress components for planar problems components   
  const MaterialProperty<Real> & _sx_alpha;
  const MaterialProperty<Real> & _sx_beta;
  
  const MaterialProperty<Real> & _sxy_alpha;
  const MaterialProperty<Real> & _sxy_beta;
  
  const MaterialProperty<Real> & _sy_alpha;
  const MaterialProperty<Real> & _sy_beta;
 
  //Input properties required are the lames constants
  const MaterialProperty<Real> &  _C11_alpha;
  const MaterialProperty<Real> &  _C11_beta;
  
  const MaterialProperty<Real> &  _C12_alpha;
  const MaterialProperty<Real> &  _C12_beta;
  
  const MaterialProperty<Real> &  _C22_alpha;
  const MaterialProperty<Real> &  _C22_beta;
  
  const MaterialProperty<Real> &  _C66_alpha;
  const MaterialProperty<Real> &  _C66_beta;
  
};
