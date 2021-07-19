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

//#ifndef BINARYCONJUGATEPHASEMATERIAL_H
//#define BINARYCONJUGATEPHASEMATERIAL_H

#pragma once
// Forward Declarations
class BinaryConjugatePhaseMaterial;

//MOOSE includes
#include "TabulatedPhaseMaterial.h"
#include "BinaryConjugatePhaseData.h"

template <>
InputParameters validParams<BinaryConjugatePhaseMaterial>();

class BinaryConjugatePhaseMaterial : public TabulatedPhaseMaterial
{
public:
  BinaryConjugatePhaseMaterial(const InputParameters & parameters);

protected:

  //Member function that returns the property value at each quadrature point
  virtual void computeQpProperties() override;
  
private:
    
    //String variable to hold the chemical potential of dep Comp A, mole fraction
    //and the inverse of the thermodyanmic factor
     std::string  _A_chem_pot_name, _xB_name, _inv_B_tf_name; //_inv_B_td_name;
    
    //Chemical potential value that this material interpolates and returns
    MaterialProperty<Real> & _A_chem_pot_val;
    
    //Mole fraction value that this material interpolates and reurns
    MaterialProperty<Real> & _xB_val;
    
    //Inv of the thermodyanamic factor that this material interpolates and returns
    MaterialProperty<Real> & _inv_B_tf_val;
    
    //Inv of the third derivative that this material interpolates and returns
    //MaterialProperty<Real> & _inv_B_td_val;
    
    //Independent variable on which this property depends
    const VariableValue & _B_diff_pot;
    
    //BinaryPhaseData is an userObject, which reads the 
    //tables for binary alloys with the properties of the phase
    //and linearly interpolates the value .
    const BinaryConjugatePhaseData & _table_object;
 
};
//#endif // BINARYCONJUGATEPHASEMATERIAL_H
