//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef MULTICOMPCHEMICALDRIVINGFORCE_H
//#define MULTICOMPCHEMICALDRIVINGFORCE_H

#pragma once

#include "ChemicalDrivingForce.h"

// Forward Declarations
class MultiCompChemicalDrivingForce;

template<>
InputParameters validParams<MultiCompChemicalDrivingForce>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class MultiCompChemicalDrivingForce : public ChemicalDrivingForce
{
public:
  MultiCompChemicalDrivingForce(const InputParameters & parameters);

protected:

  virtual Real computeValue() override;
  
  //Component C in two phases
  const VariableValue & _xC_alpha;
  const VariableValue & _xC_beta;
  
  //Diffusion potential of component C in each phase
  const MaterialProperty<Real> & _C_diff_pot_alpha;
  const MaterialProperty<Real> & _C_diff_pot_beta;
};
//#endif //MULTICOMPCHEMICALDRIVINGFORCE_H
