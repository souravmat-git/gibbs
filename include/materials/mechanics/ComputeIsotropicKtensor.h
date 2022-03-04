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
class ComputeIsotropicKtensor;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

//template <>
//InputParameters validParams<ComputeIsotropicKtensor>();

class ComputeIsotropicKtensor : public Material
{
public:
  ComputeIsotropicKtensor(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the jump in strain or the strain difference
  //and its derivative wrt to overall strain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //unit normal
    const MaterialProperty<RealVectorValue> & _n;

    //Stiffness tensor of alpha and beta phase
    const MaterialProperty<RankFourTensor> & _stiffness_alpha;
    const MaterialProperty<RankFourTensor> & _stiffness_beta;
    const MaterialProperty<Real> & _h;

    MaterialProperty<RankTwoTensor> & _K_val;
    MaterialProperty<RankTwoTensor> & _dK_dh_val;
};
