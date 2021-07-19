//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef CHEMICALDRIVINGFORCEKKS_H
//#define CHEMICALDRIVINGFORCEKKS_H
#pragma once

#include "AuxKernel.h"

// Forward Declarations
class ChemicalDrivingForceKKS;

template<>
InputParameters validParams<ChemicalDrivingForceKKS>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class ChemicalDrivingForceKKS : public AuxKernel
{
public:
  ChemicalDrivingForceKKS(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  //The gradient of a coupled variable
  const VariableValue & _diff_pot;
  
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _omega_alpha;
  const MaterialProperty<Real> & _omega_beta;
};
//#endif //CHEMICALDRIVINGFORCEKKS_H
