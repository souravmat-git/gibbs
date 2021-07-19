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
class StrainDependentTaylorApproximation;

//Include dependencies
#include "TabulatedPhaseMaterial.h"

template <>
InputParameters validParams<StrainDependentTaylorApproximation>();

 
class StrainDependentTaylorApproximation : public TabulatedPhaseMaterial
{
public:
  StrainDependentTaylorApproximation(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
  const std::string _base_name;

  const Real & _xB_eqm;
  const Real & _B_tf_eqm;
  const Real & _B_diff_pot_eqm;
  const Real & _A_chem_pot_eqm;
  
  //Derivative of eigenstrain with diffusion potential
  const MaterialProperty<RankTwoTensor>  & _stress;
  const MaterialProperty<RankTwoTensor>  & _deigen_dmu;
  const MaterialProperty<RankTwoTensor>  & _d2eigen_dmu2;
  const MaterialProperty<RankFourTensor> & _stiffness;
  
  //Compute the following properties
  MaterialProperty<Real> & _xB_val;
  MaterialProperty<Real> & _inv_B_tf_val;
  MaterialProperty<Real> & _A_chem_pot_val;
  
  const VariableValue & _B_diff_pot;  
};
