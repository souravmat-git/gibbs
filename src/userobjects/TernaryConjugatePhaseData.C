//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryConjugatePhaseData.h"
#include "BilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", TernaryConjugatePhaseData);

template <>
InputParameters
validParams<TernaryConjugatePhaseData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryConjugatePhaseData::TernaryConjugatePhaseData(const InputParameters & parameters)
  : ThermoChemicalProperties(parameters),
   _table_name(getParam<FileName>("table_name")),
   _table_reader(_table_name, &_communicator),
   _col_names(_table_reader.getNames())
{
  //Lines begining with # are comments
  _table_reader.setComment("#");  
  //tab is a delimiter
  _table_reader.setDelimiter(" ");
}

//This is a destructor
TernaryConjugatePhaseData::~TernaryConjugatePhaseData() 
{}

void
TernaryConjugatePhaseData::initialSetup()
{
  // Check to see if _file_name supplied exists. If it does, that data
  // will be used. If it does not exist, data will be generated and then
  // written to _file_name.
  std::ifstream file(_table_name.c_str());
  if (file.good())
  {
    _console << "Reading tabulated properties from " << _table_name << "\n";
    //read the data in the table and convert into double
    _table_reader.read(); 
     
    //Check the first column is always mole_fraction of components B & C
      if (_col_names[0] != "B_diff_pot" && _col_names[1] != "C_diff_pot") 
        mooseError("First and second column must be mole fraction !!"); 
   }       
   
    //Vector variable to store the array of diffusion potentials of components B &C
     _B_diff_pot = _table_reader.getData(_col_names[0]);
     
     _C_diff_pot = _table_reader.getData(_col_names[1]);
    
    //Make the data unique sort-> unique -> resize for comp B
     std::sort(_B_diff_pot.begin(), _B_diff_pot.end());
    _it_B = std::unique(_B_diff_pot.begin(),_B_diff_pot.end());
    _B_diff_pot.resize(std::distance(_B_diff_pot.begin(),_it_B));
    
    //Size of B vector which is assumed to be #of rows
    _num_xB = _B_diff_pot.size();
    
    //Make the data unique sort-> unique-> resize for comp C
     std::sort(_C_diff_pot.begin(), _C_diff_pot.end());
    _it_C = std::unique(_C_diff_pot.begin(), _C_diff_pot.end());
    _C_diff_pot.resize(std::distance(_C_diff_pot.begin(),_it_C));
    
    //Size of C vector which is assumed to #of colns
    _num_xC = _C_diff_pot.size();
      
    //Variable to store the array of chemical potential of A
    std::vector<Real> _chem_pot_A = _table_reader.getData(_col_names[2]);
    
    //Variable to store the array of mole fraction  of component B
    std::vector<Real> _xB = _table_reader.getData(_col_names[3]);
    
    //Variable to store the array of mole fraction  of component C
    std::vector<Real> _xC = _table_reader.getData(_col_names[4]);
    
    //variable to store the inverse of second derivates with respect to B
    std::vector<Real> _inv_tfactor_B = _table_reader.getData(_col_names[5]);
    
    //variable to store the inverse of second derivates with respect to BC
    std::vector<Real> _inv_tfactor_BC = _table_reader.getData(_col_names[6]);
    
    //variable to store the inverse of second derivates with respect to C
    std::vector<Real> _inv_tfactor_C = _table_reader.getData(_col_names[7]);
    
    //Declare the size of the matrix, Check PieceWiseBilinearMaterial.C
    
    _chemAMatrix.reshape(_num_xC,_num_xB);
    _xBMatrix.reshape(_num_xC, _num_xB);
    _xCMatrix.reshape(_num_xC, _num_xB);
    _inv_tfactorBMatrix.reshape(_num_xC,_num_xB);
    _inv_tfactorBCMatrix.reshape(_num_xC, _num_xB);
    _inv_tfactorCMatrix.reshape(_num_xC, _num_xB);
    
    //In case of MooseError "Reference outside of ColumnMajorMatrix bounds"
    //check the size of the matrix 
    //std::cout << _num_xC << _num_xB;
    //Assign the vector elements to the matrix elements
    for(unsigned int j=0; j < _num_xC; j++)
      for(unsigned int i=0; i < _num_xB; i++)
      {
        _chemAMatrix(i,j) = _chem_pot_A[i+_num_xC*j];
        _xBMatrix(i,j)    = _xB[i + _num_xC*j];
        _xCMatrix(i,j)    = _xC[i + _num_xC*j];
        _inv_tfactorBMatrix(i,j) = _inv_tfactor_B[i + _num_xC*j];
        _inv_tfactorBCMatrix(i,j)= _inv_tfactor_BC[i + _num_xC*j];
        _inv_tfactorCMatrix(i,j) = _inv_tfactor_C[i + _num_xC*j];
      }

    //Set the data for interpolating chemical_pot_A
   _interpolate_chem_pot_A = 
            libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot, _chemAMatrix);
    
    //Set the data for interpolating mole fraction of B
   _interpolate_xB = 
            libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot ,_xBMatrix);
    
     //Set the data for interpolating mole fraction of C
    _interpolate_xC = 
            libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot ,_xCMatrix);

    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_B = 
            libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot , _inv_tfactorBMatrix);    
       
    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_BC = 
            libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot , _inv_tfactorBCMatrix); 
    
    //Set the data for interpolating second derivatives
    _interpolate_inv_tf_C  = 
          libmesh_make_unique<BilinearInterpolation>(_B_diff_pot,_C_diff_pot , _inv_tfactorCMatrix); 
}

Real
TernaryConjugatePhaseData::A_chem_pot(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //Return the chemical potential of dependent component
  return (_interpolate_chem_pot_A->sample(_B_diff_pot,_C_diff_pot));
}

Real
TernaryConjugatePhaseData::xB(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //return 0.0;
  //Return the mole fraction of compB
  return (_interpolate_xB->sample(_B_diff_pot,_C_diff_pot));
}

Real
TernaryConjugatePhaseData::xC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //return 0.0;
  //Return the mole fraction of compB
  return (_interpolate_xC->sample(_B_diff_pot,_C_diff_pot));
}

Real
TernaryConjugatePhaseData::inv_therm_factor_B(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //return the inverse of second derivative;
  return (_interpolate_inv_tf_B->sample(_B_diff_pot, _C_diff_pot));
}

Real
TernaryConjugatePhaseData::inv_therm_factor_BC(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //return 2.0;
  return (_interpolate_inv_tf_BC->sample(_B_diff_pot, _C_diff_pot));
}

Real
TernaryConjugatePhaseData::inv_therm_factor_C(const Real& _B_diff_pot, const Real& _C_diff_pot) const
{
  //return 2.0;
  return (_interpolate_inv_tf_C->sample(_B_diff_pot, _C_diff_pot));
}
