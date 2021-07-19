//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ACBulk.h"

class GBinaryMultiPhaseDrivingForce;

template<>
InputParameters validParams<GBinaryMultiPhaseDrivingForce>();

/**
  *This class enforces the following 
  *Equation in the KKS model
  *df/dphi = dh/dphi * (omega_2 - omega_1)
  *omega_2- Grandpotential of the 2 phase
  *omega_1 - GrandPotential of the 1 phase
  *Assuming mobility is not constant
  **/

class GBinaryMultiPhaseDrivingForce :  public ACBulk<Real>
{
  public: 
    GBinaryMultiPhaseDrivingForce(const InputParameters & parameters);
  
  protected:
  
    virtual Real computeDFDOP(PFFunctionType type);
    virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  
  //One can extend these to N-phases....by adding more phase
    const VariableValue& _phase_2;
    unsigned int _phase_2_var;
    
    const VariableValue& _phase_3;
    unsigned int _phase_3_var;
  
    const VariableValue& _phase_4;
    unsigned int _phase_4_var;
  
    const VariableValue& _phase_5;
    unsigned int _phase_5_var; 

    //Coupled variable diffusion potential of comp B
    const VariableValue & _B_diff_pot;
    unsigned int _B_diff_pot_var;
    
    //string variable to hold the material property name from the input file 
    std::string  _A_chem_pot_1_name, _A_chem_pot_2_name; 
  
    //string variable to hold the material property name from the input file
    std::string  _xB_1_name, _xB_2_name;
  
    //interpolation name d2hsigma_dphitheta2
    std::string _dh_name, _d2h_name, _d2h_2_name, _d2h_3_name,_d2h_4_name, _d2h_5_name;
         
    //chemical potential of the dependent component in both phases
    const MaterialProperty<Real> & _A_chem_pot_1;
    const MaterialProperty<Real> & _A_chem_pot_2;
    
    //mu(A).(mu(B) - mu(A)) = x_B (2-1)
    const MaterialProperty<Real> & _xB_1;
    const MaterialProperty<Real> & _xB_2; 
    
    //In addition this kernel requires the first and second derivatives
    // of the interpolation function
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
  
   //Only this part of the code needs to be modified for more phases  
   //For each phase field variable in the system this will increase
    const MaterialProperty<Real> & _d2h_2;
    const MaterialProperty<Real> & _d2h_3;
    const MaterialProperty<Real> & _d2h_4;
    const MaterialProperty<Real> & _d2h_5;
    
    //A non-dimensional parameter
    const MaterialProperty<Real> & _nd_factor;   
};
