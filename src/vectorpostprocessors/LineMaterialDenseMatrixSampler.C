//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "LineMaterialDenseMatrixSampler.h"

registerMooseObject("gibbsApp", LineMaterialDenseMatrixSampler);
defineLegacyParams(LineMaterialDenseMatrixSampler);

InputParameters
LineMaterialDenseMatrixSampler::validParams()
{
  InputParameters params = LineMaterialSamplerBase<DenseMatrix<Real>>::validParams();
  params.addClassDescription("Samples real-valued material properties for all quadrature points in "
                             "all elements that are intersected by a specified line");
  params.addRequiredParam<unsigned int>("i", "The ith row of matrix C");
  params.addRequiredParam<unsigned int>("j", "The jth column of matrix C");
                             
  return params;
}

LineMaterialDenseMatrixSampler::LineMaterialDenseMatrixSampler(const InputParameters & parameters)
  : LineMaterialSamplerBase<DenseMatrix<Real>>(parameters),
  _i(getParam<unsigned int>("i")),
  _j(getParam<unsigned int>("j"))
{
}

    
Real
LineMaterialDenseMatrixSampler::getScalarFromProperty(const DenseMatrix<Real> & property, const Point & /*curr_point*/){ 
   return property(_i,_j);
}



