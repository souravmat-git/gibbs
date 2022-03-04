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

//#ifndef BINARYPHASEMATERIAL_H
//#define BINARYPHASEMATERIAL_H

#pragma once

// Forward Declarations
class BinaryPhaseMaterial;

//MOOSE includes
#include "TabulatedPhaseMaterial.h"
#include "BinaryPhaseData.h"

//template <>
//InputParameters validParams<BinaryPhaseMaterial>();

class BinaryPhaseMaterial : public TabulatedPhaseMaterial
{
public:
  BinaryPhaseMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string _f_phase_name, _B_diff_pot_phase_name, _B_therm_factor_phase_name;

    //free energy value that this material interpolates and returns
    MaterialProperty<Real> & _f_phase_val;

    //diffusion potential value that this material interpolates and reurns
    MaterialProperty<Real> & _B_diff_pot_phase_val;

     //Thermodynamic factor that this material interpolates and returns
    MaterialProperty<Real> & _B_therm_factor_phase_val;

    //Independent variable on which this property depends
    const VariableValue & _mol_fraction_B;

    //BinaryPhaseData is an userObject, which reads the
    //tables for binary alloys with the properties of the phase
    //and linearly interpolates the value .
    const BinaryPhaseData & _table_object;

};
//#endif // BINARYPHASEMATERIAL_H
