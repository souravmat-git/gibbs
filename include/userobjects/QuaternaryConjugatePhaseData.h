//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef QUATERNARYCONJUGATEPHASEDATA_H
//#define QUATERNARYCONJUGATEPHASEDATA_H

#pragma once
// Forward Declarations
class QuaternaryConjugatePhaseData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "TrilinearInterpolation.h"
#include "DelimitedFileReader.h"

template <>
InputParameters validParams<QuaternaryConjugatePhaseData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties classs
//This derived class functions is implemented in a
//similar way as PieceWise Bilinear Material.
//data of chemical potential, and mole fraction
//given by the user and then supply it to the material class
//TabulatedConjugatePhaseMaterial.

class QuaternaryConjugatePhaseData : public ThermoChemicalProperties
{
public:
  QuaternaryConjugatePhaseData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~QuaternaryConjugatePhaseData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given the diffusion potentials this function returns the chemical pot of A
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated chemical potential based on the diffusion potential of B and C
  virtual Real A_chem_pot(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;

  //Return interpolated phase composition of comp B based on the diffusion potentials
  virtual Real xB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated phase composition of comp C based on the diffusion potential
  virtual Real xC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated phase composition of comp C based on the diffusion potential
  virtual Real xD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;

  //Return interpolated inv second derivative w.r.t B
  virtual Real inv_therm_factor_B(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real inv_therm_factor_C(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated second derivative w.r.t D
  virtual Real inv_therm_factor_D(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated inv second derivative w.r.t B,C
  virtual Real inv_therm_factor_BC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated inv second derivative w.r.t B,D
  virtual Real inv_therm_factor_BD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return interpolated inv second derivative w.r.t C,D
  virtual Real inv_therm_factor_CD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  

private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<TrilinearInterpolation>_interpolate_chem_pot_A, 
                                           _interpolate_xB,
                                           _interpolate_xC,
                                           _interpolate_xD,
                                           _interpolate_inv_tf_B,
                                           _interpolate_inv_tf_C,
                                           _interpolate_inv_tf_D,
                                           _interpolate_inv_tf_BC,
                                           _interpolate_inv_tf_BD,
                                           _interpolate_inv_tf_CD;

   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    std::vector<Real> _B_diff_pot, _C_diff_pot, _D_diff_pot;
    
    //Moose utility csv_reader defined for reading the table
    //Note that DelimitedFileReader is a class defined within 
    // a namespace MoseUtils and hence to create object
    //<namespace>::<class> <object>
    MooseUtils::DelimitedFileReader _table_reader;
     
    //This reference variable holds the column names of the file
    const std::vector<std::string>& _col_names;
    
    //Declare iterators for making the vector unique
    std::vector<Real>::iterator _it_B, _it_C, _it_D;
    
    //Declare variables to store the number of data points
    unsigned int _num_xB, _num_xC, _num_xD;    
};
//#endif // QUATERNARYCONJUGATEPHASEDATA_H
