//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef TERNARYCHEMPOTENTIALDATA_H
//#define TERNARYCHEMPOTENTIALDATA_H

#pragma once
// Forward Declarations
class TernaryChemPotentialData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "BilinearInterpolation.h"
#include "DelimitedFileReader.h"
template <typename>
class ColumnMajorMatrixTempl;
typedef ColumnMajorMatrixTempl<Real> ColumnMajorMatrix;

template <>
InputParameters validParams<TernaryChemPotentialData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class TernaryChemPotentialData : public ThermoChemicalProperties
{
public:
  TernaryChemPotentialData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~TernaryChemPotentialData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section

  //Return interpolated chemical potential of component A
  virtual Real A_chem_pot(const Real& _xB, const Real& _xC) const;
 
  //Return interpolated second derivative mu(A).x(B)
  virtual Real thermodynamic_factor_AB(const Real& _xB, const Real& _xC) const;
  
  //Return interpolated second derivative mu(A).x(C)
  virtual Real thermodynamic_factor_AC(const Real& _xB, const Real& _xC) const;
 
private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<BilinearInterpolation> _interpolate_chem_pot_A,
                                           _interpolate_therm_factor_AB,
                                           _interpolate_therm_factor_AC;

   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    std::vector<Real> _xB, _xC;
    
    ColumnMajorMatrix  _chemAMatrix, _tfactorABMatrix, _tfactorACMatrix;
    
    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within 
    // a namespace MoseUtils and hence to create object
    //<namespace>::<class> <object>
    MooseUtils::DelimitedFileReader _table_reader;
     
    //This reference variable holds the column names of the file
    const std::vector<std::string>& _col_names;
    
    //Declare iterators for making the vector unique
    std::vector<Real>::iterator _it_B,_it_C;
    
    //Declare variables to store the number of data points
    MooseIndex(_xB) _num_xB;
    MooseIndex(_xC) _num_xC;    
};
//#endif // TERNARYCHEMPOTENTIALDATA_H
