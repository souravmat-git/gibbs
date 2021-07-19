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
class PrerequisiteTensors_V1;

//MOOSe includes
#include "VariableUnitNormalBase.h"

template <>
InputParameters validParams<PrerequisiteTensors_V1>();

class PrerequisiteTensors_V1 : public VariableUnitNormalBase
{
public:
  PrerequisiteTensors_V1(const InputParameters & parameters);

  //This class calculates the two required second rank tensors
  // inv_D only limited to alpha-beta regions
  //and S only limited to beta-gamma regions
  //To save the computation time, we limit the calculation to
  //the alpha-beta and beta-gamma interfacial regions

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //const MaterialProperty<RealVectorValue> & _n1; //alpha + beta
    //const MaterialProperty<RealVectorValue> & _n2; //beta + gamma

    //Stiffness tensor of alpha, beta and gamma phases
    const MaterialProperty<RankFourTensor> & _alpha_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _beta_elasticity_tensor;
    const MaterialProperty<RankFourTensor> & _gamma_elasticity_tensor;

    //interpolation functions
    const MaterialProperty<Real> & _h_alpha;
    const MaterialProperty<Real> & _h_beta;
    const MaterialProperty<Real> & _h_gamma;

    MaterialProperty<RankTwoTensor>  & _lambda_hash;
    //MaterialProperty<RankTwoTensor>  & _L_star;

    //These tensors depend on n and hence there
    //inverse blows up in the bulk
    MaterialProperty<RankTwoTensor>  & _inv_D;
    MaterialProperty<RankTwoTensor>  & _S;

};
