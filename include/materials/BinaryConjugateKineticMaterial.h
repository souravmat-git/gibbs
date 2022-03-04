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
class BinaryConjugateKineticMaterial;

//MOOSE includes
#include "TabulatedKineticMaterial.h"
#include "BinaryConjugateMobilityData.h"

//template <>
//InputParameters validParams<BinaryConjugateKineticMaterial>();

class BinaryConjugateKineticMaterial : public TabulatedKineticMaterial
{
public:
  BinaryConjugateKineticMaterial(const InputParameters & parameters);

  static InputParameters validParams();

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;

private:

    //String variable to hold L_BB and dL_BB_muB
     std::string  _L_BB_name, _dL_BB_muB_name;

    //L_BB value that this material interpolates and returns
    MaterialProperty<Real> & _L_BB_val;

    //dL_BB_muB value that this material interpolates and reurns
    MaterialProperty<Real> & _dL_BB_muB_val;

    //Independent variable on which this property depends
    const VariableValue & _B_diff_pot;

    //BinaryPhaseData is an userObject, which reads the
    //tables for binary alloys with the properties of the phase
    //and linearly interpolates the value .
    const BinaryConjugateMobilityData & _table_object;

};
