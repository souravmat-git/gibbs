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
class DerivativePrerequisiteTensors_V1;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<DerivativePrerequisiteTensors_V1>();

class DerivativePrerequisiteTensors_V1 : public Material
{
public:
  DerivativePrerequisiteTensors_V1(const InputParameters & parameters);

  //This class returns the derivatives of the D tensor
  // with respect to phase-field variables phi_alpha, phi_beta and phi_gamma

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    const MaterialProperty<RealVectorValue> & _n1; //alpha + beta
    const MaterialProperty<RealVectorValue> & _n2; //beta + gamma

    //Stiffness tensor of alpha, beta and gamma phases
    const MaterialProperty<RankFourTensor> & _alpha_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _beta_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _gamma_elasticity_tensor;

    //Derivative of interpolation fuctions
    const MaterialProperty<Real> & _dhalpha_dphialpha;
    const MaterialProperty<Real> & _dhalpha_dphibeta;
    const MaterialProperty<Real> & _dhalpha_dphigamma;

    const MaterialProperty<Real> & _dhgamma_dphialpha;
    const MaterialProperty<Real> & _dhgamma_dphibeta;
    const MaterialProperty<Real> & _dhgamma_dphigamma;

    //These tensors are needed
    const MaterialProperty<RankTwoTensor>  & _lambda_hash;
    //const MaterialProperty<RankTwoTensor>  & _lambda_star;
    //const MaterialProperty<RankTwoTensor>  & _L_hash;
    const MaterialProperty<RankTwoTensor>  & _S;

    //These tensors need to be calculated
    MaterialProperty<RankTwoTensor>  & _dD_dphialpha;
    MaterialProperty<RankTwoTensor>  & _dD_dphibeta;
    MaterialProperty<RankTwoTensor>  & _dD_dphigamma;
};
