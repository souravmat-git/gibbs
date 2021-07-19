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
class QElasticDrivingForce;

template<>
InputParameters validParams<QElasticDrivingForce>();

/**
  *Auxiliary kernel responsible for computing 
  *driving force due to solid-solid phase trannsformation
  * this assumes the driving force to be the difference 
  * between the elastic energy density between phases
  */

class QElasticDrivingForce : public AuxKernel
{
public:
  QElasticDrivingForce(const InputParameters & parameters);

protected:
  /**
    *AuxKerenls MUST override computeValue. computevalue() is called on
    *every quadrature point. For Nodal Auxiliary variables those quadrature
    *points coincide with nodes
    */
  virtual Real computeValue() override;
  
  const MaterialProperty<Real> & _h;
  const MaterialProperty<Real> & _dh;
  
  const MaterialProperty<Real> & _fel_alpha;
  const MaterialProperty<Real> & _fel_beta;
  
  const MaterialProperty<Real> & _sx_alpha;
  const MaterialProperty<Real> & _sx_beta;
  
  const MaterialProperty<Real> & _a;
  
  const MaterialProperty<Real> & _nd_factor;
};
