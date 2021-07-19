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
class AnisotropicStrainJump3P_V1;

//MOOSe includes
#include "VariableUnitNormalBase.h"

template <>
InputParameters validParams<AnisotropicStrainJump3P_V1>();

class AnisotropicStrainJump3P_V1 : public VariableUnitNormalBase
{
public:
  AnisotropicStrainJump3P_V1(const InputParameters & parameters);

  //This class calculates the inverse of K tensor

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //const MaterialProperty<RealVectorValue> & _n1;
    //const MaterialProperty<RealVectorValue> & _n2;

    //Stiffness tensor of alpha, beta and gamma phases
    const MaterialProperty<RankFourTensor> & _alpha_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _beta_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _gamma_elasticity_tensor;

    //interpolation function
    const MaterialProperty<Real> & _h_alpha;
    const MaterialProperty<Real> & _h_beta;
    const MaterialProperty<Real> & _h_gamma;

    //Fetch the prerequiste second rank tensors
    //const MaterialProperty<RankTwoTensor> & _lambda_star;
    //const MaterialProperty<RankTwoTensor> & _L_hash;
    const MaterialProperty<RankTwoTensor> & _inv_D;
    const MaterialProperty<RankTwoTensor> & _S;

    //Fetch the prerequiste vectors
    const MaterialProperty<RealVectorValue> & _m_alpha1;
    const MaterialProperty<RealVectorValue> & _m_beta1;
    const MaterialProperty<RealVectorValue> & _psi_beta2;
    const MaterialProperty<RealVectorValue> & _psi_gamma2;

    //calculate the following two vectors
    MaterialProperty<RealVectorValue>  & _a_alpha_beta;
    MaterialProperty<RealVectorValue>  & _a_beta_gamma;

};
