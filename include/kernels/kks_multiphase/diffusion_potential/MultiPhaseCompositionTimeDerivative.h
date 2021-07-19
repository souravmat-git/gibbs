//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "CoupledTimeDerivative.h"

// Forward Declaration
class MultiPhaseCompositionTimeDerivative;

template <>
InputParameters validParams<MultiPhaseCompositionTimeDerivative>();

/**
 * This calculates a modified coupled time derivative that multiplies the time derivative of a
 * coupled variable by a function of the variables
 */
class MultiPhaseCompositionTimeDerivative : public CoupledTimeDerivative
{
public:
   MultiPhaseCompositionTimeDerivative(const InputParameters & parameters);
   
protected:

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
private:

  //One can extend these to N-phases....by adding more phase
  const VariableValue & _phase_2;
  unsigned int _phase_2_var;
    
  const VariableValue & _phase_3;
  unsigned int _phase_3_var;
  
  //string variable to hold the material property name from the input file
  std::string  _xB_1_name, _xB_2_name;
  
  //string variable to hold the material property name from the input file 
  std::string  _inv_B_tf_1_name, _inv_B_tf_2_name; 
  
  
  //interpolation name 
  std::string _dh_name, _d2h_name, _d2h_2_name, _d2h_3_name;
    
  //mu(A).(mu(B) - mu(A)) = x_B (2-1)
  const MaterialProperty<Real> & _xB_1;
  const MaterialProperty<Real> & _xB_2; 
  
  //First derivative with respect to diffusion potential
  const MaterialProperty<Real> & _inv_B_tf_1;
  const MaterialProperty<Real> & _inv_B_tf_2;
    
  //In addition this kernel requires the first and second derivatives
  // of the interpolation function
  const MaterialProperty<Real> & _dh;
  const MaterialProperty<Real> & _d2h;
  
  //Only this part of the code needs to be modified for more phases  
  //For each phase field variable in the system this will increase
  const MaterialProperty<Real> & _d2h_2;
  const MaterialProperty<Real> & _d2h_3;
    
};
