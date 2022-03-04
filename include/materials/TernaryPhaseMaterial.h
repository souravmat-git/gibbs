//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//*This code is based on the code
//*from TabulatedFluidProperties.C

//#ifndef TERNARYPHASEMATERIAL_H
//#define TERNARYPHASEMATERIAL_H

#pragma once

// Forward Declarations
class TernaryPhaseMaterial;

//MOOSE includes
#include "TabulatedPhaseMaterial.h"
#include "TernaryPhaseData.h"

//template <>
//InputParameters validParams<TernaryPhaseMaterial>();

class TernaryPhaseMaterial : public TabulatedPhaseMaterial
{
public:
  TernaryPhaseMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _f_energy_name, _B_diff_pot_name, _C_diff_pot_name,
                _B_therm_factor_name, _BC_therm_factor_name,_C_therm_factor_name;

    //free energy value that this material interpolates and returns
   MaterialProperty<Real> & _f_energy_val;

    //diffusion potential of comp B that this material interpolates and reurns
   MaterialProperty<Real> & _B_diff_pot_val;

   //diffusion potential of comp C that this material interpolates and returns
   MaterialProperty<Real> & _C_diff_pot_val;

   //Thermodynamic factor w.r.t Bthat this material interpolates and returns
    MaterialProperty<Real> & _B_therm_factor_val;

    //Thermodynamic factor w.r.t BC that this material interpolates and returns
    MaterialProperty<Real> & _BC_therm_factor_val;

    //Thermodynamic factor w.r.t C that this material interpolates and returns
    MaterialProperty<Real> & _C_therm_factor_val;

    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_B;

    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_C;

    //TernaryPhaseData is an userObject, which reads the
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
    const TernaryPhaseData & _table_object;

};
//#endif // TERNARYPHASEMATERIAL_H
