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
class DerivativePrerequisiteTensors;

//MOOSe includes
#include "VariableUnitNormalBase.h"

//template <>
//InputParameters validParams<DerivativePrerequisiteTensors>();

class DerivativePrerequisiteTensors : public VariableUnitNormalBase
{
public:
  DerivativePrerequisiteTensors(const InputParameters & parameters);

  static InputParameters validParams();

  //This class returns the derivatives of the D and S tensor
  // with respect to phase-field variables phi_alpha, phi_beta and phi_gamma

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Stiffness tensors of alpha, beta and gamma phases
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

    //These tensors need to be calculated
    MaterialProperty<RankTwoTensor>  & _dD_dphialpha;
    MaterialProperty<RankTwoTensor>  & _dD_dphibeta;
    MaterialProperty<RankTwoTensor>  & _dD_dphigamma;

    MaterialProperty<RankTwoTensor> & _dS_dphialpha;
    MaterialProperty<RankTwoTensor> & _dS_dphibeta;
    MaterialProperty<RankTwoTensor> & _dS_dphigamma;
};
