//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef DIFFUSIONFLUX_H
//#define DIFFUSIONFLUX_H
#pragma once

#include "AuxKernel.h"

// Forward Declarations
class DiffusionFlux;

template<>
InputParameters validParams<DiffusionFlux>();

/**
  *Auxiliary kernel responsible for computing the flux
  *thus the concentration gradient
  */

class DiffusionFlux : public AuxKernel
{
public:
  DiffusionFlux(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

// Will hold 0,1,2 correspoding to x,y,z.
  int _component;

//The gradient of a coupled variable
  const VariableGradient  & _conc_gradient;
  Real _diffusivity;
};

//#endif //DIFFUSIONFLUX_H 

    

