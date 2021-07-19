//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef CHEMICALDRIVINGFORCE_H
//#define CHEMICALDRIVINGFORCE_H

#pragma once
#include "AuxKernel.h"

// Forward Declarations
class ChemicalDrivingForce;

template<>
InputParameters validParams<ChemicalDrivingForce>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class ChemicalDrivingForce : public AuxKernel
{
public:
  ChemicalDrivingForce(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  //Component B in two phases
  const VariableValue & _xB_alpha;
  const VariableValue & _xB_beta;
  
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _f_alpha;
  const MaterialProperty<Real> & _f_beta;
  
  //Diffusion potential of component B in each phase
  const MaterialProperty<Real> & _B_diff_pot_alpha;
  const MaterialProperty<Real> & _B_diff_pot_beta;
  
  const MaterialProperty<Real> & _nd_factor;
  
};
//#endif //CHEMICALDRIVINGFORCE_H
