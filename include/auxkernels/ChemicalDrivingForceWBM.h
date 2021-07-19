//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef CHEMICALDRIVINGFORCEWBM_H
//#define CHEMICALDRIVINGFORCEWBM_H

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class ChemicalDrivingForceWBM;

template<>
InputParameters validParams<ChemicalDrivingForceWBM>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class ChemicalDrivingForceWBM : public AuxKernel
{
public:
  ChemicalDrivingForceWBM(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;

  //Composition is a coupled variable
  const VariableValue & _comp;
  
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _free_energy_alpha;
  const MaterialProperty<Real> & _free_energy_beta;
  const MaterialProperty<Real> & _dfalpha_dc;
  const MaterialProperty<Real> & _dfbeta_dc;
};
//#endif //CHEMICALDRIVINGFORCEWBM_H
