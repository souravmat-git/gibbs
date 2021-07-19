//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This file was written by S.Chatterjee

#ifndef DRIVINGFORCEWBM_H
#define DRIVINGFORCEWBM_H

#include "Kernel.h"

class DrivingForceWBM;

template<>
InputParameters validParams<DrivingForceWBM>();

/**
  *This class enforces the following 
  *Equation in the WBM model
  *h^{\prime}(f_beta - f_alpha)= 0
  **/

class DrivingForceWBM :  public Kernel
{
  public: 
    DrivingForceWBM(const InputParameters & parameters);
  
  protected:

   //Note that the following member functions are defined 
   //virtual in the base_class Kernel.C it is not mandatory 
   //to declare it virtual. The keyword override means 
   //that the function definition is different from the base_class
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  private:
    
    const VariableValue & _comp;
    unsigned int _comp_var;
    
    //const VariableValue & _diff_pot;
    //unsigned int _diff_pot_var;
    
    //Keep in mind the free energy are functions of 
    //composition and not phase compositions
    //Note the difference in the naming of the two variables 
    //dfalpha_dc and df_dcalpha
    const MaterialProperty<Real> & _dh;
    const MaterialProperty<Real> & _d2h;
    const MaterialProperty<Real> & _free_energy_alpha;
    const MaterialProperty<Real> & _dfalpha_dc;
    const MaterialProperty<Real> & _free_energy_beta;
    const MaterialProperty<Real> & _dfbeta_dc;
    const MaterialProperty<Real> & _L;
    const MaterialProperty<Real> & _nd_factor;
};
#endif // DIFFUSIONPOTENTIALWBM_H
