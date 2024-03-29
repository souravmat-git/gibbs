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
class ComputeSmallStrainAlpha;

//MOOSe includes
#include "ComputeStrainBase.h"


template <>
InputParameters validParams<ComputeSmallStrainAlpha>();

class ComputeSmallStrainAlpha : public ComputeStrainBase
{
public:
  ComputeSmallStrainAlpha(const InputParameters & parameters);

  //This class calculates the strain for the alpha phase


protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Jump in (compatible) strain depends on the homogenization
    const MaterialProperty<RankTwoTensor>  & _strain_jump;

    //Interpolation function which depends on phase-field
    const MaterialProperty<Real> & _h;

};
