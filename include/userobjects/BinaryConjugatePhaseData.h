//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef BINARYCONJUGATEPHASEDATA_H
//#define BINARYCONJUGATEPHASEDATA_H

#pragma once
// Forward Declarations
class BinaryConjugatePhaseData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "LinearInterpolation.h"
#include "DelimitedFileReader.h"

template <>
InputParameters validParams<BinaryConjugatePhaseData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class BinaryConjugatePhaseData : public ThermoChemicalProperties
{
public:
  BinaryConjugatePhaseData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~BinaryConjugatePhaseData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated chemical potential of the dependent comp
  virtual Real A_chem_pot(const Real& _B_diff_pot) const;

  //Return interpolated diffusion potential of component B
  virtual Real xB(const Real& _B_diff_pot) const;

  //Return interpolated second derivative
  virtual Real inv_therm_factor_B(const Real& _B_diff_pot) const;
  
  //Return interpolated third derivative
  //virtual Real inv_third_deriv_B(const Real& _B_diff_pot) const;

private:
   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    //variable to hold interpolated chemical potential of A
    std::unique_ptr<LinearInterpolation> _interpolate_chem_pot_A;
    
    //variable to hold interpolated chem_pot of comp A And B
    std::unique_ptr<LinearInterpolation> _interpolate_xB;
    
    //variable to hold interpolated inv  of thermodynamic factor;
    std::unique_ptr<LinearInterpolation> _interpolate_inv_tf_B;
    
    //variable to hold interpolated inv of the third derivative
    //std::unique_ptr<LinearInterpolation> _interpolate_inv_td_B;
       
    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within 
    // a namespace MoseUtils and hence to create object
    //<namespace>::<class> <object>
    MooseUtils::DelimitedFileReader _table_reader;
     
    //This reference variable holds the column names of the file
    const std::vector<std::string>& _col_names;
    
};
//#endif // BINARYCONJUGATEPHASEDATA_H
