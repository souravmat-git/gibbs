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
class ComputeQDrivingForceOffDiag;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//template <>
//InputParameters validParams<ComputeQDrivingForceOffDiag>();

class ComputeQDrivingForceOffDiag : public Material
{
public:
  ComputeQDrivingForceOffDiag(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    const MaterialProperty<RankTwoTensor> & _alpha_stress;
    const MaterialProperty<RankTwoTensor> & _beta_stress;

    //Stiffness tensor of alpha and beta phase
    const MaterialProperty<RankFourTensor> & _alpha_stiffness;
    const MaterialProperty<RankFourTensor> & _beta_stiffness;

    //Strain jump and its derivatives
    const MaterialProperty<RankTwoTensor>  & _strain_jump;
    const MaterialProperty<RankTwoTensor>  & _dstrainjump_dphi;
    const MaterialProperty<RankFourTensor> & _ds_de;

    //Derivative interpolation
    const MaterialProperty<Real> & _h;
    const MaterialProperty<Real> & _dh;

    //Non-dimensional factor
    const MaterialProperty<Real> & _nd_factor;

    MaterialProperty<RankTwoTensor>  & _d2Fdcdstrain;
};
