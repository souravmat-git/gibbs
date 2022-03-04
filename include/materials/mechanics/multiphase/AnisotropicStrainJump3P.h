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
class AnisotropicStrainJump3P;

//MOOSe includes
#include "VariableUnitNormalBase.h"

//template <>
//InputParameters validParams<AnisotropicStrainJump3P>();

class AnisotropicStrainJump3P : public VariableUnitNormalBase
{
public:
  AnisotropicStrainJump3P(const InputParameters & parameters);

  static InputParameters validParams();

  //This class returns the jump vectors
  //a_alpha_beta and a_beta_gamma

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Fetch the prerequiste second rank tensors
    const MaterialProperty<RankTwoTensor> & _inv_D;
    const MaterialProperty<RankTwoTensor> & _S;

    //Fetch the prerequiste vectors in alpha-beta regions
    const MaterialProperty<RealVectorValue> & _m_alpha1;
    const MaterialProperty<RealVectorValue> & _m_beta1;
    const MaterialProperty<RealVectorValue> & _Z1;

    //Fetch the prerequsite vectors in the beta-gamma regions
    const MaterialProperty<RealVectorValue> & _psi_beta2;
    const MaterialProperty<RealVectorValue> & _psi_gamma2;
    const MaterialProperty<RealVectorValue> & _Z2;

    //calculate the following two vectors
    MaterialProperty<RealVectorValue>  & _a_alpha_beta;
    MaterialProperty<RealVectorValue>  & _a_beta_gamma;

};
