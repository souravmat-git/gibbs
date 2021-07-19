//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef TOTALENERGY_H
//#define TOTALENERGY_H

#pragma once
#include "InterfaceEnergy.h"

// Forward Declarations
class TotalEnergy;

template<>
InputParameters validParams<TotalEnergy>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class TotalEnergy : public InterfaceEnergy
{
public:
 TotalEnergy(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;
  

  const MaterialProperty<Real> & _f_alpha;
  const MaterialProperty<Real> & _f_beta;
  //const MaterialProperty<Real> & _h_alpha;
  const MaterialProperty<Real> & _h_beta;
 
};
//#endif //TOTALENERGY_H
