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
class AnisotropicStrainJumpMaterialN;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankThreeTensor.h"
#include "RankFourTensor.h"

//template <>
//InputParameters validParams<AnisotropicStrainJumpMaterialN>();

class AnisotropicStrainJumpMaterialN : public Material
{
public:
  
  AnisotropicStrainJumpMaterialN(const InputParameters & parameters);
  static InputParameters validParams();

  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

  //unit normal
  //const MaterialProperty<RealVectorValue> & _n;
  const VariableGradient & _grad_eta;

  //Number of displacements
  unsigned int _ndisp;

  //displacement variables
  //std::vector<const VariableValue *> _disp;

  //Gradients of displacements as a vector
  std::vector<const VariableGradient *> _grad_disp;

  //Eigenstrains of alpha and beta phase
  const MaterialProperty<RankTwoTensor> & _eigen_alpha;
  const MaterialProperty<RankTwoTensor> & _eigen_beta;

  //Stiffness tensor of alpha and beta phase
  const MaterialProperty<RankFourTensor> & _stiffness_alpha;
  const MaterialProperty<RankFourTensor> & _stiffness_beta;

  //inverse of K tensor
  const MaterialProperty<RankTwoTensor>  & _inv_K;


  MaterialProperty<RealVectorValue> & _a_val;
};
