//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef BINARYCONJUGATEMOBILITYDATA_H
//#define BINARYCONJUGATEMOBILITYDATA_H

#pragma once
// Forward Declarations
class BinaryConjugateMobilityData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "LinearInterpolation.h"
#include "DelimitedFileReader.h"


template <>
InputParameters validParams<BinaryConjugateMobilityData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class BinaryConjugateMobilityData : public ThermoChemicalProperties
{
public:
  BinaryConjugateMobilityData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~BinaryConjugateMobilityData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated mobility as a function of diffusion potential
  virtual Real L_BB(const Real& _B_diff_pot) const;

  //Return the first derivative of the mobility with respect to diff pot
  virtual Real dL_BB_muB(const Real& _B_diff_pot) const;

private:
   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    //variable to hold interpolated mobility L_BB
    std::unique_ptr<LinearInterpolation> _interpolate_L_BB;
    
    //variable to hold interpolated derivative
    std::unique_ptr<LinearInterpolation> _interpolate_dL_BB_muB;
       
    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within 
    // a namespace MoseUtils and hence to create object
    //<namespace>::<class> <object>
    MooseUtils::DelimitedFileReader _table_reader;
     
    //This reference variable holds the column names of the file
    const std::vector<std::string>& _col_names;
    
};
//#endif // BINARYCONJUGATEMOBILITYDATA_H
