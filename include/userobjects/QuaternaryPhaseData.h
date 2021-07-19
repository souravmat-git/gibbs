//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef QUATERNARYPHASEDATA_H
#define QUATERNARYPHASEDATA_H

// Forward Declarations
class QuaternaryPhaseData;

//MOOSE includes
#include "ThermoChemicalProperties.h"
#include "TrilinearInterpolation.h"
#include "DelimitedFileReader.h"


template <>
InputParameters validParams<QuaternaryPhaseData>();

//TabulatedPhaseData is a derived class from the
//ThermoChemicalProperties class s
//This derived class functions is implemented in a
// similar way as TabulatedFluidProperties.
//data of free_energy, diffusion_potential and second derivative
//given by the user and then supply it to the material class
//TabulatedPhaseMaterial.

class QuaternaryPhaseData : public ThermoChemicalProperties
{
public:
  QuaternaryPhaseData(const InputParameters & parameters);

  //Destructor is a member function which has the
  //same name as the class but they do not take any argument
  //nor return anything source: WSatvich
  virtual ~QuaternaryPhaseData(); //Destructor

  virtual void initialSetup() override;

  //Given a mole fraction this function returns the free energy
  //Note that the const is placed at the end of a
  //member function declartion !! Ex:BlockAverageValue.h
  //These member functions are called by the object in the Material section

  //Return interpolated diffusion potential of component B
  virtual Real B_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated diffusion potential of component C
  virtual Real C_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated diffusion potential of component C
  virtual Real D_diff_pot(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated second derivative mu(B).x(B)
  virtual Real thermodynamic_factor_B(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated second derivative mu(C).x(C)
  virtual Real thermodynamic_factor_C(const Real& _xB, const Real& _xC, const Real& _xD) const;

    //Return interpolated second derivative mu(D).x(D)
  virtual Real thermodynamic_factor_D(const Real& _xB, const Real& _xC, const Real& _xD) const;

    //Return interpolated second derivative mu(B).x(C)
  virtual Real thermodynamic_factor_BC(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated second derivative mu(B).x(D)
  virtual Real thermodynamic_factor_BD(const Real& _xB, const Real& _xC, const Real& _xD) const;

  //Return interpolated second derivative mu(C).x(D)
  virtual Real thermodynamic_factor_CD(const Real& _xB, const Real& _xC, const Real& _xD) const;

private:

    //Code based on PiecewiseBilinear
    std::unique_ptr<TrilinearInterpolation>_interpolate_chem_pot_A,
                                           _interpolate_chem_pot_B,
                                           _interpolate_chem_pot_C,
                                           _interpolate_chem_pot_D,
                                           _interpolate_therm_factor_B,
                                           _interpolate_therm_factor_C,
                                           _interpolate_therm_factor_D,
                                           _interpolate_therm_factor_BC,
                                           _interpolate_therm_factor_BD,
                                           _interpolate_therm_factor_CD;


    //variable to hold the string type table name from input file
    FileName _table_name;

    //Three independent variables for a quartenary system
    std::vector<Real> _xB, _xC, _xD;


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
    unsigned int _num_xB, _num_xD; //_num_xC   
};
#endif // QUATERNARYPHASEDATA_H
