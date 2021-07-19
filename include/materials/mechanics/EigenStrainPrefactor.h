//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
class EigenStrainPrefactor;

//Moose inludes
#include "Material.h"

template <>
InputParameters validParams<EigenStrainPrefactor>();

class EigenStrainPrefactor : public Material
{
public:
  EigenStrainPrefactor(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  
private:
     
   //base_name prepended to material name
   const std::string _base_name;     
   
   //Gas constant and temperature
   //required for non-dimensionalizaton
   const Real & _R;
   const Real & _T; 
      
   //Parameters required to calculate the diffusion potential
   const Real & _xB_eqm;
   const Real & _B_tf_eqm;
   const Real & _B_diff_pot_eqm;
      
   const Real & _xB_o;
     
   MaterialProperty<Real> & _prefactor;
   MaterialProperty<Real> & _dprefactor;
   MaterialProperty<Real> & _d2prefactor;
   
   const VariableValue & _B_diff_pot;

};
