//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef TERNARYCONJUGATEPHASEDATA_H
//#define TERNARYCONJUGATEPHASEDATA_H

#pragma once
// Forward Declarations
class TernaryConjugatePhaseData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "BilinearInterpolation.h"
#include "DelimitedFileReader.h"
template <typename>
class ColumnMajorMatrixTempl;
typedef ColumnMajorMatrixTempl<Real> ColumnMajorMatrix;

template <>
InputParameters validParams<TernaryConjugatePhaseData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties classs
//This derived class functions is implemented in a
//similar way as PieceWise Bilinear Material.
//data of chemical potential, and mole fraction
//given by the user and then supply it to the material class
//TabulatedConjugatePhaseMaterial.

class TernaryConjugatePhaseData : public ThermoChemicalProperties
{
public:
  TernaryConjugatePhaseData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~TernaryConjugatePhaseData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated chemical potential based on the diffusion potential of B and C
  virtual Real A_chem_pot(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

  //Return interpolated phase composition of comp B based on the diffusion potentials
  virtual Real xB(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated phase composition of comp C based on the diffusion potential
  virtual Real xC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

  //Return interpolated inv second derivative w.r.t B
  virtual Real inv_therm_factor_B(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated inv second derivative w.r.t B,C
  virtual Real inv_therm_factor_BC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real inv_therm_factor_C(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<BilinearInterpolation> _interpolate_chem_pot_A, 
                                           _interpolate_xB,
                                           _interpolate_xC,
                                           _interpolate_inv_tf_B,
                                           _interpolate_inv_tf_BC,
                                           _interpolate_inv_tf_C;

   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    std::vector<Real> _B_diff_pot, _C_diff_pot;
    
    ColumnMajorMatrix _chemAMatrix, _xBMatrix, _xCMatrix, 
                      _inv_tfactorBMatrix, _inv_tfactorBCMatrix, _inv_tfactorCMatrix;
    
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
    unsigned int _num_xB, _num_xC;    
};
//#endif // TERNARYCONJUGATEPHASEDATA_H
