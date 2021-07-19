//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MatrixFunctionIC.h"

// MOOSE includes
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "DelimitedFileReader.h"

//C++ includes
#include <cmath>
#include <fstream>

registerMooseObject("gibbsApp", MatrixFunctionIC);

defineLegacyParams(MatrixFunctionIC);


InputParameters
MatrixFunctionIC::validParams()
{
  InputParameters params = Function::validParams();
  params.addRequiredParam<FileName>("table_name", "Table with x y r");
  params.addRequiredParam<unsigned int>("nx", "Number of finite element in X-direction");
  params.addRequiredParam<unsigned int>("ny", "Number of finite element in Y-directionr");
  return params;
}

MatrixFunctionIC::MatrixFunctionIC(const InputParameters & parameters)
  : Function(parameters),
   _fe_problem(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
   _mesh(_fe_problem.mesh()),
   //obtain this from the file
   _nx(getParam<unsigned int>("nx")),
   _ny(getParam<unsigned int>("ny")),
   _table_name(getParam<FileName>("table_name")),
   _table_reader(_table_name, &_communicator),
   _data(_table_reader.getData())
{
  //Check if number of elements supplied are correct
  unsigned int _nelem  = _nx *_ny;
  if(_nelem != _mesh.nElem())
    mooseError("Total mesh elements != Mesh elements in generated mesh");

  //Check if file is okay
  std::ifstream file(_table_name.c_str());
  if (!file.good())
    mooseError("Can't read the input file format");
  else
    _table_reader.read();

  //Lines begining with # are comments
  _table_reader.setComment("#");
  //comma is a delimiter
  _table_reader.setDelimiter(",");

  //Check the size of data and the number of nodes
  unsigned int _nrows = _data[0].size();
  unsigned int _ncols = _data.size();

  unsigned int _datasize  = _nrows *_ncols;

   if (_datasize != _mesh.nNodes())
    mooseWarning("total number of nodes != data size \n"
                  "The variable value will be simply set to zero");
}

Real
MatrixFunctionIC::value(Real /*t*/, const Point & pt) const
{

  const Real x = pt(0);
  const Real y = pt(1);

  const Real xmin = _mesh.getMinInDimension(0);
  const Real xmax = _mesh.getMaxInDimension(0);

  const Real ymin  = _mesh.getMinInDimension(1);
  const Real ymax  = _mesh.getMaxInDimension(1);

  const Real dx = (xmax - xmin)/_nx;
  const Real dy = (ymax - ymin)/_ny;

  unsigned int i = std::round(x/dx);
  unsigned int j = std::round(y/dy);

  if (i < _data[0].size() && j < _data.size())
    return _data[i][j];
  else
    return 0;

}
