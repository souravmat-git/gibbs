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
class StrainJump_3P;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

template <>
InputParameters validParams<StrainJump_3P>();

class StrainJump_3P : public Material
{
public:
  StrainJump_3P(const InputParameters & parameters);

  //This class calculates the inverse of K tensor

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //const VariableGradient & _gradphi_alpha;
    //const VariableGradient & _gradphi_beta ;
    //const VariableGradient & _gradphi_gamma;

    //Gradients of displacements as a vector
    //std::vector<const VariableGradient *> _grad_disp;
    const MaterialProperty<RankTwoTensor> & _total_strain;

    //Stiffness tensor of alpha, beta and gamma phases
    const MaterialProperty<Real> & _mat_const_alpha;
    const MaterialProperty<Real> & _mat_const_beta;
    const MaterialProperty<Real> & _mat_const_gamma;

    //interpolation function
    const MaterialProperty<Real> & _h_alpha;
    const MaterialProperty<Real> & _h_beta;
    const MaterialProperty<Real> & _h_gamma;

    MaterialProperty<Real>  & _a_alpha_beta;
    MaterialProperty<Real>  & _a_beta_gamma;

};
