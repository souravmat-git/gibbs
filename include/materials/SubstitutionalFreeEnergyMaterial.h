//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef SUBSTITUTIONALFREEENERGYMATERIAL_H
//#define SUBSTITUTIONALFREEENERGYMATERIAL_H

#pragma once

#include "Material.h"
#include "cmath"

//Forward declaration
class SubstitutionalFreeEnergyMaterial;

template <>
InputParameters validParams<SubstitutionalFreeEnergyMaterial>();

/**
  * This class defines
  * the substitutional free energy material
  * it takes two input parameters 
  * GHSERNi and GHSERAl
  */

class SubstitutionalFreeEnergyMaterial : public Material
{
  public:
    SubstitutionalFreeEnergyMaterial(const InputParameters & parameters);

  protected:
    virtual void computeQpProperties() override;
    //can be accessed by base members, and derived class 
       
    /*Material constants are protected member variables*
     and thus can be used in derived classes*/
    const Real _Ec;
    //characteristic energy (Units: J/m^3)
    const Real _Vm;
    //Molar volume of phase (Units: mol/m^3)
    const Real _R;
    // gas constant (Units: J/molK)
    const Real & _T;
    // temperatue of the simulation (Units: K)
    const Real & _GHSERA;
    // GHSER of component A at 298K and 1atm pressure (Units: J)
    const Real & _GHSERB;
    // GHSER of component B at 298K & 1atm pressure
    const std::vector<std::string> _RK_names;
    const std::vector<Real> _RK_val;
    //Redlich-Kister expansion terms

   private:
      //component B is assumed to be independent
     const VariableValue & _xB;
      
     std::string _free_energy_name, _B_diff_pot_name, _chi_name;

     MaterialProperty<Real> & _free_energy;
     //Returns the free energy of a binary A-B alloy

     MaterialProperty<Real> & _B_diff_pot;
     //Returns the diffusion potential of component B

     MaterialProperty<Real> & _chi;
    //Returns the second derivative with resoect to B
    //In the literature, this term is called the thermodynamic factor
 
     Real _G_ref, _G_conf,_G_ex; //computing the free energy
     Real _dG_ref, _dG_conf, _dG_ex1, _dG_ex2;
     Real _d2G_ref, _d2G_conf, _d2G_ex1, _d2G_ex2, _d2G_ex3;
};
//#endif //SUBSTITUTIONALFREEENERGYMATERIAL_H
