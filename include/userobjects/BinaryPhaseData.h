//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef BINARYPHASEDATA_H
//#define BINARYPHASEDATA_H

#pragma once
// Forward Declarations
class BinaryPhaseData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "LinearInterpolation.h"
#include "DelimitedFileReader.h"


template <>
InputParameters validParams<BinaryPhaseData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class BinaryPhaseData : public ThermoChemicalProperties
{
public:
  BinaryPhaseData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~BinaryPhaseData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated free energy
  virtual Real free_energy(const Real& _xB) const;

  //Return interpolated diffusion potential of component B
  virtual Real B_diff_pot(const Real& _xB) const;

  //Return interpolated second derivative
  virtual Real thermodynamic_factor(const Real& _xB) const;

private:
   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    //variable to hold interpolated free energy
    std::unique_ptr<LinearInterpolation> _interpolate_free_energy;
    
    //variable to hold interpolated chem_pot of comp A And B
    std::unique_ptr<LinearInterpolation> _interpolate_B_diff_pot;
    
    //variable to hold interpolated thermodynamic factor;
    std::unique_ptr<LinearInterpolation> _interpolate_B_therm_factor;
       
    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within 
    // a namespace MoseUtils and hence to create object
    //<namespace>::<class> <object>
    MooseUtils::DelimitedFileReader _table_reader;
     
    //This reference variable holds the column names of the file
    const std::vector<std::string>& _col_names;
    
};
//#endif // BINARYPHASEDATA_H
