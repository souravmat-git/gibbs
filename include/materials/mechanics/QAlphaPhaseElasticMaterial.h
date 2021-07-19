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
class QAlphaPhaseElasticMaterial;

//MOOSe includes
#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "MooseMesh.h"

template <>
InputParameters validParams<QAlphaPhaseElasticMaterial>();

class QAlphaPhaseElasticMaterial : public Material
{
public:
  QAlphaPhaseElasticMaterial(const InputParameters & parameters);
  
  //This class calculates the strain, stress and the elastic energy
  //of a phase with or without a prescribed eigenstrain

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:

    //********************** Input ***************************//  
    
    //First get the name of the phase
    const std::string _phase_name;  
      
    //Total strain
    const MaterialProperty<RankTwoTensor> & _total_strain;
     
    //Eigenstrain or the transformation strain   
    const MaterialProperty<RankTwoTensor>  & _alpha_eigen;
    
    //Jump in (compatible) strain depends on the homogenization
    const MaterialProperty<RankTwoTensor>  & _strain_jump;
    
    //Stiffness tensor
    const MaterialProperty<RankFourTensor> & _alpha_stiffness;
    
    //Interpolation function which depends on phase-field
    const MaterialProperty<Real> & _h;
    
    //********************** Output ***************************//
       
    //Phase elastic strain
    MaterialProperty<RankTwoTensor> & _elastic_strain_val;
    
    //Stress in the x-direction
    MaterialProperty<RankTwoTensor> & _stress_val;
    
    //Jacobian
    //MaterialProperty<RankFourTensor> & _Jacobian_mult;
    
    //Elastic strain energy
    MaterialProperty<Real> & _elastic_energy_val;    
};
