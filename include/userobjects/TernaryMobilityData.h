//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#ifndef TERNARYMOBILITYDATA_H
//#define TERNARYMOBILITYDATA_H

#pragma once
// Forward Declarations
class TernaryMobilityData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "BilinearInterpolation.h"
#include "DelimitedFileReader.h"

template <typename>
class ColumnMajorMatrixTempl;
typedef ColumnMajorMatrixTempl<Real> ColumnMajorMatrix;

template <>
InputParameters validParams<TernaryMobilityData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class TernaryMobilityData : public ThermoChemicalProperties
{
public:
  TernaryMobilityData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~TernaryMobilityData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated free energy depdenent of comp B and comp C
  virtual Real L_BB(const Real& _xB, const Real& _xC) const;

  //Return interpolated diffusion potential of component B
  virtual Real L_BC(const Real& _xB, const Real& _xC) const;
  
  //Return interpolated diffusion potential of component C
  virtual Real L_CC(const Real& _xB, const Real& _xC) const;

  //Return interpolated second derivative w.r.t B
  virtual Real dL_BB_xB(const Real& _xB, const Real& _xC) const;
  
  //Return interpolated second derivative w.r.t B,C
  virtual Real dL_BB_xC(const Real& _xB, const Real& _xC) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real dL_BC_xB(const Real & _xB, const Real & _xC) const;
  
  virtual Real dL_BC_xC(const Real & _xB, const Real & _xC) const;
  
  //Return interpolated second derivative w.r.t C
  virtual Real dL_CC_xB(const Real & _xB, const Real & _xC) const;
  
  virtual Real dL_CC_xC(const Real & _xB, const Real & _xC) const;

private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<BilinearInterpolation> _interpolate_LBB, 
                                           _interpolate_LBC,
                                           _interpolate_LCC,
                                           _interpolate_dLBB_xB,
                                           _interpolate_dLBB_xC,
                                           _interpolate_dLBC_xB,
                                           _interpolate_dLBC_xC,
                                           _interpolate_dLCC_xB,
                                           _interpolate_dLCC_xC;
                                           
   
    //variable to hold the string type table name from input file
    FileName _table_name; 
    
    std::vector<Real> _xB, _xC;
    
    ColumnMajorMatrix _LBBMatrix,    _LBCMatrix,    _LCCMatrix,
                      _dLBBxBMatrix, _dLBBxCMatrix, _dLBCxBMatrix,
                      _dLBCxCMatrix, _dLCCxBMatrix, _dLCCxCMatrix;
    
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
//#endif // TERNARYMOBILITYDATA_H
