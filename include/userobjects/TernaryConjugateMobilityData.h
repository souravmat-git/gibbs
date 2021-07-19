//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef TERNARYCONJUGATEMOBILITYDATA_H
//#define TERNARYCONJUGATEMOBILITYDATA_H

#pragma once
// Forward Declarations
class TernaryConjugateMobilityData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "BilinearInterpolation.h"
#include "DelimitedFileReader.h"

template <typename>
class ColumnMajorMatrixTempl;
typedef ColumnMajorMatrixTempl<Real> ColumnMajorMatrix;

template <>
InputParameters validParams<TernaryConjugateMobilityData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class TernaryConjugateMobilityData : public ThermoChemicalProperties
{
public:
  TernaryConjugateMobilityData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~TernaryConjugateMobilityData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated free energy depdenent of comp B and comp C
  virtual Real L_BB(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

  //Return interpolated diffusion potential of component B
  virtual Real L_BC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated diffusion potential of component C
  virtual Real L_CC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

  //Return interpolated second derivative w.r.t B
  virtual Real dL_BB_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated second derivative w.r.t B,C
  virtual Real dL_BB_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real dL_BC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  virtual Real dL_BC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real dL_CC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot) const;
  
  virtual Real dL_CC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot) const;

private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<BilinearInterpolation> _interpolate_LBB, 
                                           _interpolate_LBC,
                                           _interpolate_LCC,
                                           _interpolate_dLBB_muB,
                                           _interpolate_dLBB_muC,
                                           _interpolate_dLBC_muB,
                                           _interpolate_dLBC_muC,
                                           _interpolate_dLCC_muB,
                                           _interpolate_dLCC_muC;
                                           
   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    std::vector<Real> _B_diff_pot, _C_diff_pot;
    
    ColumnMajorMatrix _LBBMatrix,    _LBCMatrix,    _LCCMatrix,
                      _dLBBmuBMatrix, _dLBBmuCMatrix, _dLBCmuBMatrix,
                      _dLBCmuCMatrix, _dLCCmuBMatrix, _dLCCmuCMatrix;
    
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
//#endif // TERNARYCONJUGATEMOBILITYDATA_H
