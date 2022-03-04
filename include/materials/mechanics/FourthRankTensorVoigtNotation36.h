//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class FourthRankTensorVoigtNotation36;

//MOOSE includes
#include "Material.h"

//template <>
//InputParameters validParams<FourthRankTensorVoigtNotation36>();

class FourthRankTensorVoigtNotation36 : public Material
{
public:
  FourthRankTensorVoigtNotation36(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    const MaterialProperty<RankFourTensor> & _C;
    //Given stress components
    std::string _C_matrix_name;
    MaterialProperty<DenseMatrix<Real>> & _V;

};
