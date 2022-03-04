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
class PlaneElasticityStrainEnergy;

//MOOSe includes
#include "Material.h"

//template <>
//InputParameters validParams<PlaneElasticityStrainEnergy>();

class PlaneElasticityStrainEnergy : public Material
{
public:
  PlaneElasticityStrainEnergy(const InputParameters & parameters);

  static InputParameters validParams();

  //This class calculates the strain from the displacement gradient
  //Using the strain-displacement relations.

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //Independent variables on which this property depends
    // are the displacement gradient vector along x
    // and the displacement gradient vector along y

    std::string _phase_name;
    std::string _fel_name;

    //Given stress and strain declared in the inpput block
    const MaterialProperty<Real> & _sxx;
    const MaterialProperty<Real> & _syy;
    const MaterialProperty<Real> & _sxy;
    const MaterialProperty<Real> & _exx;
    const MaterialProperty<Real> & _eyy;
    const MaterialProperty<Real> & _exy;

    MaterialProperty<Real> & _fel_val;
};
