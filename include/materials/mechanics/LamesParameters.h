//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// Forward Declarations
class LamesParameters;

//MOOSe includes
#include "Material.h"
#include "RankFourTensor.h"

//template <>
//InputParameters validParams<LamesParameters>();

class LamesParameters : public Material
{
public:
  LamesParameters(const InputParameters & parameters);

  static InputParameters validParams();
  //This class returns the lames parameter from th stiffness tensor

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

  const std::string _phase_name;

  //Stiffness tensor of alpha and beta phase
  const MaterialProperty<RankFourTensor> & _stiffness;

  MaterialProperty<Real> & _lambda;
  MaterialProperty<Real> & _mu;

};
