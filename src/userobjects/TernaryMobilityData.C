//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TernaryMobilityData.h"
#include "BilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", TernaryMobilityData);

template <>
InputParameters
validParams<TernaryMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

TernaryMobilityData::TernaryMobilityData(const InputParameters & parameters)
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
TernaryMobilityData::~TernaryMobilityData() 
{}

void
TernaryMobilityData::initialSetup()
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
      if (_col_names[0] != "x_B" && _col_names[1] != "x_C") 
        mooseError("First and second column must be mole fraction !!"); 
   }       
   
    //Vector variable to store the array of mole_fraction of components B &C
     _xB = _table_reader.getData(_col_names[0]);
     
     _xC = _table_reader.getData(_col_names[1]);
    
    //Make the data unique sort-> unique -> resize for comp B
     std::sort(_xB.begin(), _xB.end());
    _it_B = std::unique(_xB.begin(), _xB.end());
    _xB.resize(std::distance(_xB.begin(),_it_B));
    
    //Size of B vector which is assumed to be #of rows
    _num_xB = _xB.size();
    
    //Make the data unique sort-> unique-> resize for comp C
     std::sort(_xC.begin(), _xC.end());
    _it_C = std::unique(_xC.begin(), _xC.end());
    _xC.resize(std::distance(_xC.begin(),_it_C));
    
    //Size of C vector which is assumed to #of colns
    _num_xC = _xC.size();
      
    //Variable to store the array of free energy
    std::vector<Real> _LBB = _table_reader.getData(_col_names[2]);
    
    //Variable to store the array of chemical potential of component A
    std::vector<Real> _LBC = _table_reader.getData(_col_names[3]);
    
    //variable to store the array of chemical potential of component B
    std::vector<Real> _LCC = _table_reader.getData(_col_names[4]);
    
    //variable to store the array of chemical potential of component C
    std::vector<Real> _dLBB_xB = _table_reader.getData(_col_names[5]);
    
    //variable to store second derivates with respect to B
    std::vector<Real> _dLBB_xC = _table_reader.getData(_col_names[6]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _dLBC_xB = _table_reader.getData(_col_names[7]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _dLBC_xC = _table_reader.getData(_col_names[8]);
    
    //variable to store second derivates with respect to BC
    std::vector<Real> _dLCC_xB = _table_reader.getData(_col_names[9]);
    
    //variable to store second derivates with respect to C
    std::vector<Real> _dLCC_xC = _table_reader.getData(_col_names[10]);
    
    
    //Declare the size of the matrix, Check PieceWiseBilinearMaterial.C
    
    _LBBMatrix.reshape(_num_xC, _num_xB);
    _LBCMatrix.reshape(_num_xC,_num_xB);
    _LCCMatrix.reshape(_num_xC, _num_xB);
    //derivative of LBB
    _dLBBxBMatrix.reshape(_num_xC, _num_xB);
    _dLBBxCMatrix.reshape(_num_xC,_num_xB);
    //derivatives of LBC
    _dLBCxBMatrix.reshape(_num_xC, _num_xB);
    _dLBCxCMatrix.reshape(_num_xC,_num_xB);
    //derivatives of LCC
    _dLCCxBMatrix.reshape(_num_xC, _num_xB);
    _dLCCxCMatrix.reshape(_num_xC,_num_xB);
    //In case of MooseError "Reference outside of ColumnMajorMatrix bounds"
    //check the size of the matrix 
    //std::cout << _num_xC << _num_xB;
    //Assign the vector elements to the matrix elements
    for(MooseIndex(_num_xC) j=0; j < _num_xC; j++)
      for(MooseIndex(_num_xB) i=0; i < _num_xB; i++)
      {
        _LBBMatrix(i,j) = _LBB[i + _num_xC*j];
        _LBCMatrix(i,j) = _LBC[i+_num_xC*j];
        _LCCMatrix(i,j) = _LCC[i + _num_xC*j];
        _dLBBxBMatrix(i,j) = _dLBB_xB[i + _num_xC*j];
        _dLBBxCMatrix(i,j) = _dLBB_xC[i + _num_xC*j];
        _dLBCxBMatrix(i,j) = _dLBC_xB[i + _num_xC*j];
        _dLBCxCMatrix(i,j) = _dLBC_xC[i + _num_xC*j];
        _dLCCxBMatrix(i,j) = _dLCC_xB[i + _num_xC*j];
        _dLCCxCMatrix(i,j) = _dLCC_xC[i + _num_xC*j];
      }
    
    //Set the data for interpolating free energy
   _interpolate_LBB = libmesh_make_unique<BilinearInterpolation>(_xB, _xC, _LBBMatrix);

    //Set the data for interpolating chemical_pot_A
   _interpolate_LBC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _LBCMatrix);
    
    //Set the data for interpolating chemical_pot_B
   _interpolate_LCC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC, _LCCMatrix);
    
    //Set the data for interpolating for L_BB derivatives with respect to comp B,C
    _interpolate_dLBB_xB = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLBBxBMatrix);

    _interpolate_dLBB_xC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLBBxCMatrix);    
     
    //Set the data for interpolating for L_BC derivatives with respect to compB, C
    _interpolate_dLBC_xB = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLBCxBMatrix);
    
    _interpolate_dLBC_xC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLBCxCMatrix);
    
    //Set the data for interpolating for L_CC derivatives with respect to comp C
    _interpolate_dLCC_xB = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLCCxBMatrix);
    
    _interpolate_dLCC_xC = libmesh_make_unique<BilinearInterpolation>(_xB,_xC,_dLCCxCMatrix);
}

Real 
TernaryMobilityData::L_BB(const Real& _xB, const Real& _xC) const
{   
    return (_interpolate_LBB->sample(_xB,_xC));   
}

Real
TernaryMobilityData::L_BC(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_LBC->sample(_xB,_xC));
}

Real
TernaryMobilityData::L_CC(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_LCC->sample(_xB,_xC));
}

Real 
TernaryMobilityData::dL_BB_xB(const Real& _xB, const Real& _xC) const
{   
    return (_interpolate_dLBB_xB->sample(_xB,_xC));   
}

Real
TernaryMobilityData::dL_BB_xC(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_dLBB_xC->sample(_xB,_xC));
}

Real
TernaryMobilityData::dL_BC_xB(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_dLBC_xB->sample(_xB,_xC));
}

Real 
TernaryMobilityData::dL_BC_xC(const Real& _xB, const Real& _xC) const
{   
    return (_interpolate_dLBC_xC->sample(_xB,_xC));   
}

Real
TernaryMobilityData::dL_CC_xB(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_dLCC_xB->sample(_xB,_xC));
}

Real
TernaryMobilityData::dL_CC_xC(const Real& _xB, const Real& _xC) const
{
  return (_interpolate_dLCC_xC->sample(_xB,_xC));
}
