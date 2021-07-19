//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryMobilityData.h"
#include "TrilinearInterpolation.h"
#include "MooseUtils.h"

registerMooseObject("gibbsApp", QuaternaryMobilityData);

template <>
InputParameters
validParams<QuaternaryMobilityData>()
{
  InputParameters params = validParams<ThermoChemicalProperties>();
  params.addRequiredParam<FileName>("table_name", "Phase data in a table");
  params.addClassDescription("Given any tabulated property for a phase..."
                              "this class returns the interpolated values"); 
  return params;
}

QuaternaryMobilityData::QuaternaryMobilityData(const InputParameters & parameters)
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
QuaternaryMobilityData::~QuaternaryMobilityData() 
{}

void
QuaternaryMobilityData::initialSetup()
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
   _xD = _table_reader.getData(_col_names[2]);
       
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

    //Size of C vector which is assumed to #of rows
    _num_xC = _xC.size();

    //Make the data unique sort-> unique-> resize for comp D
     std::sort(_xD.begin(), _xD.end());
    _it_D = std::unique(_xD.begin(), _xD.end());
    _xD.resize(std::distance(_xD.begin(),_it_D));

    //Size of D vector which is assumed to #of rows
    _num_xD = _xD.size();
        

    //Variable to store the array of diagonal coeffecients of the mobility L
    std::vector<Real> _LBB = _table_reader.getData(_col_names[3]);
    std::vector<Real> _LCC = _table_reader.getData(_col_names[4]);
    std::vector<Real> _LDD = _table_reader.getData(_col_names[5]);

    //Variable to store the array of off-diagonal coeffecients of the mobility L
    std::vector<Real> _LBC = _table_reader.getData(_col_names[6]);
    std::vector<Real> _LBD = _table_reader.getData(_col_names[7]);
    std::vector<Real> _LCD = _table_reader.getData(_col_names[8]);
    
    //variable to store the first derivative of each coeff. with respect to B
    std::vector<Real> _dLBB_xB = _table_reader.getData(_col_names[9]);
    std::vector<Real> _dLCC_xB = _table_reader.getData(_col_names[10]);
    std::vector<Real> _dLDD_xB = _table_reader.getData(_col_names[11]);
   
    std::vector<Real> _dLBC_xB = _table_reader.getData(_col_names[12]);
    std::vector<Real> _dLBD_xB = _table_reader.getData(_col_names[13]);
    std::vector<Real> _dLCD_xB = _table_reader.getData(_col_names[14]);
    
    //variable to store the first derivative of each coeff. with respect to C
    std::vector<Real> _dLBB_xC = _table_reader.getData(_col_names[15]);
    std::vector<Real> _dLCC_xC = _table_reader.getData(_col_names[16]);
    std::vector<Real> _dLDD_xC = _table_reader.getData(_col_names[17]);
   
    std::vector<Real> _dLBC_xC = _table_reader.getData(_col_names[18]);
    std::vector<Real> _dLBD_xC = _table_reader.getData(_col_names[19]);
    std::vector<Real> _dLCD_xC = _table_reader.getData(_col_names[20]);
    
    //variable to store the first derivative of each coeff. with respect to D
    std::vector<Real> _dLBB_xD = _table_reader.getData(_col_names[21]);
    std::vector<Real> _dLCC_xD = _table_reader.getData(_col_names[22]);
    std::vector<Real> _dLDD_xD = _table_reader.getData(_col_names[23]);
   
    std::vector<Real> _dLBC_xD = _table_reader.getData(_col_names[24]);
    std::vector<Real> _dLBD_xD = _table_reader.getData(_col_names[25]);
    std::vector<Real> _dLCD_xD = _table_reader.getData(_col_names[26]);
        
     //Set the data for interpolating diagonal terms of the mobility matrix
   _interpolate_LBB = 
        libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD, _LBB); 
   _interpolate_LCC = 
        libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD, _LCC);     
   _interpolate_LDD = 
        libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD, _LDD);
         
    //Set the data for interpolating off-diagonal terms of the mobility matrix
   _interpolate_LBC = 
        libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD, _LBC);
   _interpolate_LBD = 
        libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD, _LBD);
   _interpolate_LCD = 
         libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD, _LCD);
         
    //Set the data for interpolating for first derivatives with respect to mole fraction of B
    _interpolate_dLBB_xB = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLBB_xB);
    _interpolate_dLCC_xB = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLCC_xB);    
    _interpolate_dLDD_xB = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_dLDD_xB);
    
    _interpolate_dLBC_xB = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBC_xB);
    _interpolate_dLBD_xB = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBD_xB);
    _interpolate_dLCD_xB = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLCD_xB);
          
     //Set the data for interpolating for first derivatives with respect to mole fraction of C
    _interpolate_dLBB_xC = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLBB_xC);
    _interpolate_dLCC_xC = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLCC_xC);    
    _interpolate_dLDD_xC = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_dLDD_xC);
    
    _interpolate_dLBC_xC = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBC_xC);
    _interpolate_dLBD_xC = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBD_xC);
    _interpolate_dLCD_xC = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLCD_xC);
          
     //Set the data for interpolating for first derivatives with respect to mole fraction of D
    _interpolate_dLBB_xD = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLBB_xD);
    _interpolate_dLCC_xD = 
         libmesh_make_unique<TrilinearInterpolation>(_xB, _xC, _xD,_dLCC_xD);    
    _interpolate_dLDD_xD = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC, _xD,_dLDD_xD);
    
    _interpolate_dLBC_xD = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBC_xD);
    _interpolate_dLBD_xD = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLBD_xD);
    _interpolate_dLCD_xD = 
          libmesh_make_unique<TrilinearInterpolation>(_xB,_xC,_xD,_dLCD_xD);   
}

//Diagonal terms
Real 
QuaternaryMobilityData::L_BB(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_LBB->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::L_CC(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_LCC->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::L_DD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_LDD->sample(_xB,_xC,_xD));
}

//Off-diagonal terms
Real 
QuaternaryMobilityData::L_BC(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_LBC->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::L_BD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_LBD->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::L_CD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_LCD->sample(_xB,_xC,_xD));
}

//First derivative with respect to xB
Real 
QuaternaryMobilityData::dL_BB_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBB_xB->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_CC_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCC_xB->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_DD_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLDD_xB->sample(_xB,_xC,_xD));
}

Real 
QuaternaryMobilityData::dL_BC_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBC_xB->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_BD_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLBD_xB->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_CD_xB(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCD_xB->sample(_xB,_xC,_xD));
}

//First derivative with respect to xC
Real 
QuaternaryMobilityData::dL_BB_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBB_xC->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_CC_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCC_xC->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_DD_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLDD_xC->sample(_xB,_xC,_xD));
}

Real 
QuaternaryMobilityData::dL_BC_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBC_xC->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_BD_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLBD_xC->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_CD_xC(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCD_xC->sample(_xB,_xC,_xD));
}

//First derivative with respect to xD
Real 
QuaternaryMobilityData::dL_BB_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBB_xD->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_CC_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCC_xD->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_DD_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLDD_xD->sample(_xB,_xC,_xD));
}

Real 
QuaternaryMobilityData::dL_BC_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{   
    return (_interpolate_dLBC_xD->sample(_xB,_xC,_xD));   
}
Real
QuaternaryMobilityData::dL_BD_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLBD_xD->sample(_xB,_xC,_xD));
}
Real
QuaternaryMobilityData::dL_CD_xD(const Real& _xB, const Real& _xC, const Real& _xD) const
{
  return (_interpolate_dLCD_xD->sample(_xB,_xC,_xD));
}
