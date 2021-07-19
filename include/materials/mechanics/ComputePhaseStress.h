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
class ComputePhaseStress;

//MOOSe includes
#include "ComputeStressBase.h"


template <>
InputParameters validParams<ComputePhaseStress>();

class ComputePhaseStress : public ComputeStressBase
{
public:
  ComputePhaseStress(const InputParameters & parameters);
  
  //This class calculates the stress and defines the elastic strain

protected:
  virtual void computeQpStress() override;
  
private:
  
    //Stiffness tensor
    const MaterialProperty<RankFourTensor> & _stiffness;
    
};
