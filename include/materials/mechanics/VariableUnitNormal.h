//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class VariableUnitNormal;
//MOOSE includes
#include "Material.h"

//template <>
//InputParameters validParams<VariableUnitNormal>();

class VariableUnitNormal : public Material
{
public:
  VariableUnitNormal(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

private:

    const VariableGradient & _grad_eta;

    //A user-defined value for the bulk phases
    const RealVectorValue _n;

    //Name of the variable
    const std::string _unit_normal_name;

    //A user-specified tolerance
    //const Real  & _tol;

    //Declare the unit normal vector
    MaterialProperty<RealVectorValue> &  _n_val;
    //and its deriavtive with respect to phase-field
    MaterialProperty<RankTwoTensor>   & _dn_dgradphi;

};
