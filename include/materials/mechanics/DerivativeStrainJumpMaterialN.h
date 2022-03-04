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
class DerivativeStrainJumpMaterialN;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

//template <>
//InputParameters validParams<DerivativeStrainJumpMaterialN>();

class DerivativeStrainJumpMaterialN : public Material
{
public:
  DerivativeStrainJumpMaterialN(const InputParameters & parameters);

  static InputParameters validParams();

   //Calculates the derivative of strain jump with respect to strain
   //and phase-field variable

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

  //unit normal
  //const MaterialProperty<RealVectorValue> & _n;
  const VariableGradient & _grad_eta;

  //a = Jump in strain
  const MaterialProperty<RealVectorValue> & _a;

  //Eigenstrains of alpha and beta phase
  const MaterialProperty<RankTwoTensor> & _eigen_alpha;
  const MaterialProperty<RankTwoTensor> & _eigen_beta;

  //Stiffness tensor of alpha and beta phase
  const MaterialProperty<RankFourTensor> & _stiffness_alpha;
  const MaterialProperty<RankFourTensor> & _stiffness_beta;

  //inverse of K and its derivative
  const MaterialProperty<RankTwoTensor>  & _inv_K;
  const MaterialProperty<RankTwoTensor>  & _dK_dh;

  //Derivative of interpolation function
  const MaterialProperty<Real> & _dh;

  MaterialProperty<RankFourTensor>  & _ds_de_val;
  MaterialProperty<RealVectorValue> & _da_dphi_val;
};
