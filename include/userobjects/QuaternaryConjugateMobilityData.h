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
class QuaternaryConjugateMobilityData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "TrilinearInterpolation.h"
#include "DelimitedFileReader.h"

template <>
InputParameters validParams<QuaternaryConjugateMobilityData>();

//TabulatedPhaseData is a derived class from the 
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class QuaternaryConjugateMobilityData : public ThermoChemicalProperties
{
public:
  QuaternaryConjugateMobilityData(const InputParameters & parameters);
  
  //Destructor is a member function which has the 
  //same name as the class but they do not take any argument
  //nor return anything source
  virtual ~QuaternaryConjugateMobilityData(); //Destructor
 
  virtual void initialSetup() override;
   
  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a 
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section
  
  //Return interpolated diagonal components of the mobility matrix
  virtual Real L_BB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real L_CC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real L_DD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return the off-diagonal components of the mobility matrix
  virtual Real L_BC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;  
  virtual Real L_BD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real L_CD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return the fisrt dervative of the interpolated mobility  w.r.t mu B
  virtual Real dL_BB_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_DD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  virtual Real dL_BC_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_BD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CD_muB(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  //Return the fisrt dervative of the interpolated mobility  w.r.t mu C
  virtual Real dL_BB_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_DD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
 
  virtual Real dL_BC_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const; 
  virtual Real dL_BD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CD_muC(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
     
  //Return the fisrt dervative of the interpolated mobility  w.r.t mu D
  virtual Real dL_BB_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CC_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_DD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
  virtual Real dL_BC_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_BD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  virtual Real dL_CD_muD(const Real& _B_diff_pot, const Real& _C_diff_pot, const Real& _D_diff_pot) const;
  
private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<TrilinearInterpolation>_interpolate_LBB, _interpolate_LCC, _interpolate_LDD,
                                           _interpolate_LBC, _interpolate_LBD, _interpolate_LCD,
                                           _interpolate_dLBB_muB, _interpolate_dLCC_muB, _interpolate_dLDD_muB,
                                           _interpolate_dLBC_muB, _interpolate_dLBD_muB, _interpolate_dLCD_muB,
                                           _interpolate_dLBB_muC, _interpolate_dLCC_muC, _interpolate_dLDD_muC,
                                           _interpolate_dLBC_muC, _interpolate_dLBD_muC, _interpolate_dLCD_muC,
                                           _interpolate_dLBB_muD, _interpolate_dLCC_muD, _interpolate_dLDD_muD,
                                           _interpolate_dLBC_muD, _interpolate_dLBD_muD, _interpolate_dLCD_muD;

                                           
   
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
    std::vector<Real>::iterator _it_B,_it_C, _it_D;
    
    //Declare variables to store the number of data points
    unsigned int _num_xB, _num_xC, _num_xD;    
};
