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
class ComputeQDrivingForceOffDiagBase_3P;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<ComputeQDrivingForceOffDiagBase_3P>();

class ComputeQDrivingForceOffDiagBase_3P : public Material
{
public:
  ComputeQDrivingForceOffDiagBase_3P(const InputParameters & parameters);

  static InputParameters validParams();

  //This is the base class to implement the Jacobian
  //corresponding to the driving force with respect to strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

  //Declare \sigma_{ij}^{\alpha}*(\frac{\partial
  //\epsilon_{ij}^{\alpha}}{\partial \epsilon_{mn}})

  RankTwoTensor dsalpha() const;
  RankTwoTensor dsbeta()  const;
  RankTwoTensor dsgamma() const;

  const MaterialProperty<RankTwoTensor> & _alpha_stress;
  const MaterialProperty<RankTwoTensor> & _beta_stress;
  const MaterialProperty<RankTwoTensor> & _gamma_stress;

  //Strain jumps
  const MaterialProperty<RankTwoTensor> & _strain_jump_alpha_beta;
  const MaterialProperty<RankTwoTensor> & _strain_jump_beta_gamma;

  //Interpolation function
  const MaterialProperty<Real> & _h_alpha;
  const MaterialProperty<Real> & _h_beta;
  const MaterialProperty<Real> & _h_gamma;

  //Strain jumps derivative with respect to strain
  const MaterialProperty<RankFourTensor> & _ds_alpha_beta_de;
  const MaterialProperty<RankFourTensor> & _ds_beta_gamma_de;

  //Overall stress
  const MaterialProperty<RankTwoTensor> & _stress;

  //Jacobian mult
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  //Non-dimensional factor
  const MaterialProperty<Real> & _nd_factor;

};
