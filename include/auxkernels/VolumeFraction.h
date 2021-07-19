//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef VOLUMEFRACTION_H
//#define VOLUMEFRACTION_H

#pragma once

#include "AuxKernel.h"
#include "ElementAverageValue.h"

// Forward Declarations
class VolumeFraction;

template<>
InputParameters validParams<ElementAverageValue>();

/**
  *Auxiliary kernel for computing 
  *the volume fraction or eta*eta
  */

class VolumeFraction : public ElementAverageValue
{
public:
  VolumeFraction(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeQpIntegral() override;
  
};
//#endif //VOLUMEFRACTION_H
