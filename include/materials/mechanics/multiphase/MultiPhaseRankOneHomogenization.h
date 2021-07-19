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
class MultiPhaseRankOneHomogenization;

//MOOSe includes
#include "VariableUnitNormalBase.h"

template <>
InputParameters validParams<MultiPhaseRankOneHomogenization>();

class MultiPhaseRankOneHomogenization : public VariableUnitNormalBase
{
public:
  MultiPhaseRankOneHomogenization(const InputParameters & parameters);

  //This class calculates the difference in compatible elastic strain
  //and supplies this value to other material class
  //Its current value depends on the local value of the jump in displacement
  //and the unit normal at that point

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Magnitude of strain jump
    const MaterialProperty<RealVectorValue> & _avec_alpha_beta;
    const MaterialProperty<RealVectorValue> & _avec_beta_gamma;

    //and their derivatives wrt phi_alpha, phi_beta, and phi_gamma
    const MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphialpha;
    const MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphibeta;
    const MaterialProperty<RealVectorValue> & _davec_alpha_beta_dphigamma;

    const MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphialpha;
    const MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphibeta;
    const MaterialProperty<RealVectorValue> & _davec_beta_gamma_dphigamma;


    //Compute two second-rank tensors re;ated to strain jump at
    //the alpha-beta and beta/gamma interface
    MaterialProperty<RankTwoTensor> & _strain_jump_alpha_beta;
    MaterialProperty<RankTwoTensor> & _strain_jump_beta_gamma;

    //Derivative of strain jumps with respect to phi_alpha
    MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphialpha;
    MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphialpha;

    //Derivative of strain jumps with respect to phi_beta
    MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphibeta;
    MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphibeta;

    //Derivative of strain jumps with respect to phi_gamma
    MaterialProperty<RankTwoTensor> & _ds_alpha_beta_dphigamma;
    MaterialProperty<RankTwoTensor> & _ds_beta_gamma_dphigamma;
};
