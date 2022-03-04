//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class ConstantMaterialModulus;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<ConstantMaterialModulus>();

class ConstantMaterialModulus : public Material
{
public:
  ConstantMaterialModulus(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the material constant for 1D elastic problems
  //given that we know E, nu.

protected:
  virtual void computeQpProperties() override;

private:

    const std::string _phase_name;

    //std::string  _phase_name;
    const MaterialProperty<RankFourTensor> & _stiffness;

    //Material constant for 1D problems
    MaterialProperty<Real> & _mat_const;
};
