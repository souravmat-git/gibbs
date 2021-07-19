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
class ElasticProperties3PGamma;

//MOOSe includes
#include "Material.h"

template <>
InputParameters validParams<ElasticProperties3PGamma>();

class ElasticProperties3PGamma : public Material
{
public:
  ElasticProperties3PGamma(const InputParameters & parameters);

  //This class calculates the strain for the alpha phase
  //for a 3-phase system

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //First get the name of the phase
    const std::string _phase_name;

    //Total strain
    const MaterialProperty<RankTwoTensor> & _total_strain;

    //Jump in (compatible) strain depends on the homogenization
    const MaterialProperty<Real> & _a_alpha_beta;
    const MaterialProperty<Real> & _a_beta_gamma;

    //Interpolation function which depends on phase-field
    const MaterialProperty<Real> & _h_alpha;
    const MaterialProperty<Real> & _h_beta;

    //Modulus of phase
    const MaterialProperty<Real> & _mat_const;

    //Phase elastic strain
    MaterialProperty<Real> & _elastic_strain_val;
    MaterialProperty<Real> & _elastic_stress_val;
    MaterialProperty<Real> & _strain_energy_val;
};
