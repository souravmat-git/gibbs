//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef BULKENERGY_H
//#define BULKENERGY_H

#pragma once
#include "AuxKernel.h"

// Forward Declarations
class BulkEnergy;

template<>
InputParameters validParams<BulkEnergy>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class BulkEnergy : public AuxKernel
{
public:
 BulkEnergy(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

   //Note h_alpha is not required since
   //h_alpha = (1-h_beta) always
  //const MaterialProperty<Real> & _h_alpha;
  const MaterialProperty<Real> & _h_beta;
  const MaterialProperty<Real> & _f_alpha;
  const MaterialProperty<Real> & _f_beta;
  
};
//#endif //BULKENERGY_H
