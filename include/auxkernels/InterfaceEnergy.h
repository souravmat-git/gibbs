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
class InterfaceEnergy;

template<>
InputParameters validParams<InterfaceEnergy>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class InterfaceEnergy : public AuxKernel
{
public:
 InterfaceEnergy(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;
  
  const VariableValue & _phase_alpha;
  const VariableGradient & _grad_phase_alpha;
  
  const VariableValue & _phase_beta;
  const VariableGradient & _grad_phase_beta;
  
  const MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _BH;
  
};
//#endif //INTERFACEENERGY_H
