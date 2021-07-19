//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef GRANDPOTENTIAL_H
//#define GRANDPOTENTIAL_H

#pragma once

#include "AuxKernel.h"
#include "ChemicalDrivingForce.h"

// Forward Declarations
class GrandPotential;

template<>
InputParameters validParams<GrandPotential>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class GrandPotential : public AuxKernel
{
public:
  GrandPotential(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  //Component B in two phases
  const VariableValue & _xB;
  
  //Component C in two phases
  const VariableValue & _xC;

  const MaterialProperty<Real> & _f;
  
  //Diffusion potential of component B in a  phase
  const MaterialProperty<Real> & _B_diff_pot;
  
  //Diffusion potential of component C in each phase
  const MaterialProperty<Real> & _C_diff_pot;
};
//#endif //GRANDPOTENTIAL_H
