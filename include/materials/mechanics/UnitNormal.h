//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class UnitNormal;
//MOOSE includes
#include "Material.h"


template <>
InputParameters validParams<UnitNormal>();

class UnitNormal : public Material
{
public:
  UnitNormal(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
private:
    
    const RealVectorValue _n;
    const std::string _unit_normal_name;
    
    //Declare the unit normal vector
    MaterialProperty<RealVectorValue> &  _n_val;
};
