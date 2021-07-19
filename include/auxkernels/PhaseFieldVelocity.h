//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


//#ifndef PHASEFIELDVELOCITY_H
//#define PHASEFIELDVELOCITY_H

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class PhaseFieldVelocity;

template<>
InputParameters validParams<PhaseFieldVelocity>();

/**
  *Auxiliary kernel responsible for computing the flux
  *thus the concentration gradient
  */

class PhaseFieldVelocity : public AuxKernel
{
public:
  PhaseFieldVelocity(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  // Will hold 0,1,2 correspoding to x,y,z.
  int _component;
  
  const VariableValue & _eta;
  
  const VariableValue & _dot_eta;
  
  const VariableGradient & _grad_eta;

};
//#endif //PHASEFIELDVELOCITY_
