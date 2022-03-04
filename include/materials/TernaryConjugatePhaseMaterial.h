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

//#ifndef TERNARYCONJUGATEPHASEMATERIAL_H
//#define TERNARYCONJUGATEPHASEMATERIAL_H

#pragma once

// Forward Declarations
class TernaryConjugatePhaseMaterial;

//MOOSE includes
#include "TabulatedPhaseMaterial.h"
#include "TernaryConjugatePhaseData.h"

//template <>
//InputParameters validParams<TernaryConjugatePhaseMaterial>();

class TernaryConjugatePhaseMaterial : public TabulatedPhaseMaterial
{
public:
  TernaryConjugatePhaseMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //String variable to hold th free energy, diffusion potential
    //and the thermodynamic facor name obtained from the input file
    std::string  _A_chem_pot_name, _xB_name, _xC_name,
                 _inv_B_tf_name,  _inv_BC_tf_name, _inv_C_tf_name;



    //chemical potential of comp A that this material interpolates and reurns
   MaterialProperty<Real> & _A_chem_pot_val;

   //Mole fraction of component B that this material interpolates and returns
   MaterialProperty<Real> & _xB_val;

   //Mole fraction of component C that this material interpolates and returns
   MaterialProperty<Real> & _xC_val;

   //Thermodynamic factor w.r.t Bthat this material interpolates and returns
    MaterialProperty<Real> & _inv_B_tf_val;

    //Thermodynamic factor w.r.t BC that this material interpolates and returns
    MaterialProperty<Real> & _inv_BC_tf_val;

    //Thermodynamic factor w.r.t C that this material interpolates and returns
    MaterialProperty<Real> & _inv_C_tf_val;

    //Independent variable on which this property depends
    const VariableValue & _B_diff_pot;

    //Independent variable on which this property depends
    const VariableValue & _C_diff_pot;

    //TernaryPhaseData is an userObject, which reads the
    //tables for ternary alloys with the properties of the phase
    //and linearly interpolates the value .
    const TernaryConjugatePhaseData & _table_object;

};
//#endif // TERNARYCONJUGATEPHASEMATERIAL_H
