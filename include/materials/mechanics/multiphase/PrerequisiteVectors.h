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
class PrerequisiteVectors;

//MOOSe includes
#include "VariableUnitNormalBase.h"

//template <>
//InputParameters validParams<PrerequisiteVectors>();

class PrerequisiteVectors : public VariableUnitNormalBase
{
public:
  PrerequisiteVectors(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Total strain
    const MaterialProperty<RankTwoTensor> & _total_strain;

    //Eigenstrains
    const MaterialProperty<RankTwoTensor> & _alpha_eigen;
    const MaterialProperty<RankTwoTensor> & _beta_eigen;
    const MaterialProperty<RankTwoTensor> & _gamma_eigen;

    //Stiffness tensor of alpha, beta and gamma phases
    const MaterialProperty<RankFourTensor> & _alpha_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _beta_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _gamma_elasticity_tensor;

    //Only needed in the alpha-beta regions
    MaterialProperty<RealVectorValue>  & _m_alpha1;
    MaterialProperty<RealVectorValue>  & _m_beta1;
    MaterialProperty<RealVectorValue>  & _Z1;

    //only needed in the beta-gamma regions
    MaterialProperty<RealVectorValue>  & _psi_beta2;
    MaterialProperty<RealVectorValue>  & _psi_gamma2;
    MaterialProperty<RealVectorValue>  & _Z2;
};
