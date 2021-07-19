//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html#pragma once
// Forward Declarations

class PartialRankOneHomogenization;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"

template <>
InputParameters validParams<PartialRankOneHomogenization>();

class PartialRankOneHomogenization : public Material
{
public:
  PartialRankOneHomogenization(const InputParameters & parameters);
  
  //This class calculates the difference in compatible elastic strain
  //and supplies this value to other material class
  //Its current value depends on the local value of the jump in displacement
  //and the unit normal at that point

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //Magnitude of strain jump
    const MaterialProperty<RealVectorValue> & _a;
    
    //and its derivative wrt h
    const MaterialProperty<RealVectorValue> & _da_dphi;
    
    //UnitNormal
    //const MaterialProperty<RealVectorValue> & _n;
    const VariableGradient & _grad_eta;

    //Compute the second-rank tensor
    MaterialProperty<RankTwoTensor> & _strain_jump;  
    MaterialProperty<RankTwoTensor> & _dstrainjump_dphi;        
};
