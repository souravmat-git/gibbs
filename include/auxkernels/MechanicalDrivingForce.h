//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class MechanicalDrivingForce;

template<>
InputParameters validParams<MechanicalDrivingForce>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-solid phase trannsformation
  * this assumes the driving force to be the difference 
  * between the elastic energy density between phases
  */

class MechanicalDrivingForce : public AuxKernel
{
public:
  MechanicalDrivingForce(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  
private:

  Real driving_force_aux() const;
  
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  const MaterialProperty<Real> & _alpha_strain_energy;
  const MaterialProperty<Real> & _beta_strain_energy;
  
  const MaterialProperty<RankTwoTensor> & _alpha_stress;
  const MaterialProperty<RankTwoTensor> & _beta_stress;
  
  const MaterialProperty<RankTwoTensor> & _strain_jump;
  
  const MaterialProperty<Real> & _nd_factor;

};
