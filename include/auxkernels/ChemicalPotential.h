//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef CHEMICALPOTENTIAL_H
//#define CHEMICALPOTENTIAL_H
#pragma once

#include "AuxKernel.h"

// Forward Declarations
class ChemicalPotential;

template<>
InputParameters validParams<ChemicalPotential>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-liquid phase trannsformation
  */

class ChemicalPotential : public AuxKernel
{
public:
  ChemicalPotential(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;
  

  const MaterialProperty<Real> & _A_chem_pot_alpha;
  const MaterialProperty<Real> & _A_chem_pot_beta;
  
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _nd_factor;
  
};
//#endif //CHEMICALPOTENTIAL_H
