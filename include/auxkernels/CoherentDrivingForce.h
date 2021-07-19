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
class CoherentDrivingForce;

template<>
InputParameters validParams<CoherentDrivingForce>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-solid phase trannsformation
  * this assumes the driving force to be the difference 
  * between the elastic energy density between phases
  */

class CoherentDrivingForce : public AuxKernel
{
public:
  CoherentDrivingForce(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  //Composition is a coupled variable
  
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _comp_energy_alpha;
  const MaterialProperty<Real> & _comp_energy_beta;
  const MaterialProperty<Real> & _A_chem_pot_alpha;
  const MaterialProperty<Real> & _A_chem_pot_beta;
};
