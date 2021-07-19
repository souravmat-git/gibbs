//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef ONSAGERFLUX_H
//#define ONSAGERFLUX_H
#pragma once

#include "AuxKernel.h"

// Forward Declarations
class OnsagerFlux;

template<>
InputParameters validParams<OnsagerFlux>();

/**
  *Auxiliary kernel responsible for computing the flux
  *thus the concentration gradient
  */

class OnsagerFlux : public AuxKernel
{
public:
  OnsagerFlux(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  // Will hold 0,1,2 correspoding to x,y,z.
  int _component;

  //The gradient of diffusion potential of comp. B
  const VariableGradient  & _B_diff_pot_gradient;
  
  //Gradient of diffusion potential of comp.C
  const VariableGradient & _C_diff_pot_gradient;
  
  const MaterialProperty<Real> & _M_11;
  const MaterialProperty<Real> & _M_12;
};
//#endif //ONSAGERFLUX_H
