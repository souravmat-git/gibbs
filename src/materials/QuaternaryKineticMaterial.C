//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "QuaternaryKineticMaterial.h"
//Only inlcude the header file of the class
//and not its dependency
registerMooseObject("gibbsApp", QuaternaryKineticMaterial);

//template <>
InputParameters
QuaternaryKineticMaterial::validParams()
{
  InputParameters params = TabulatedKineticMaterial::validParams();
  params.addRequiredParam<MaterialPropertyName>("L_BB", "Onsager mobility L_BB");
  params.addRequiredParam<MaterialPropertyName>("L_CC", "Onsager mobility L_CC");
  params.addRequiredParam<MaterialPropertyName>("L_DD", "Onsager mobility L_DD");
  params.addRequiredParam<MaterialPropertyName>("L_BC", "Onsager mobility L_BC");
  params.addRequiredParam<MaterialPropertyName>("L_BD", "Onsager mobility L_BD");
  params.addRequiredParam<MaterialPropertyName>("L_CD", "Onsager mobility L_CD");
  //w.r.t B
  params.addParam<MaterialPropertyName>("dL_BB_xB",0.0,"L_BB w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CC_xB",0.0,"L_CC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_DD_xB",0.0,"L_DD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BC_xB",0.0,"L_BC w.r.t B");
  params.addParam<MaterialPropertyName>("dL_BD_xB",0.0,"L_BD w.r.t B");
  params.addParam<MaterialPropertyName>("dL_CD_xB",0.0,"L_CD w.r.t B");
  //w.r.t C
  params.addParam<MaterialPropertyName>("dL_BB_xC",0.0,"L_BB w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CC_xC",0.0,"L_CC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_DD_xC",0.0,"L_DD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BC_xC",0.0,"L_BC w.r.t C");
  params.addParam<MaterialPropertyName>("dL_BD_xC",0.0,"L_BD w.r.t C");
  params.addParam<MaterialPropertyName>("dL_CD_xC",0.0,"L_CD w.r.t C");
  //w.r.t D
  params.addParam<MaterialPropertyName>("dL_BB_xD",0.0,"L_BB w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CC_xD",0.0,"L_CC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_DD_xD",0.0,"L_DD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BC_xD",0.0,"L_BC w.r.t D");
  params.addParam<MaterialPropertyName>("dL_BD_xD",0.0,"L_BD w.r.t D");
  params.addParam<MaterialPropertyName>("dL_CD_xD",0.0,"L_CD w.r.t D");
  //Coupled variables
  params.addCoupledVar("xB", "Mole fraction of component B");
  params.addCoupledVar("xC", "Mole fraction of component C");
  params.addCoupledVar("xD", "Mole fraction of component D");
  params.addRequiredParam<UserObjectName>("table_object", "Name of the table with phase properties");
  params.addClassDescription("Given any tabulated property data for a phase..."
                              "this class returns the interpolated values");
  return params;
}

QuaternaryKineticMaterial::QuaternaryKineticMaterial(const InputParameters & parameters)
  : TabulatedKineticMaterial(parameters),
   _L_BB_name(getParam<MaterialPropertyName>("L_BB")),
   _L_CC_name(getParam<MaterialPropertyName>("L_CC")),
   _L_DD_name(getParam<MaterialPropertyName>("L_DD")),
   _L_BC_name(getParam<MaterialPropertyName>("L_BC")),
   _L_BD_name(getParam<MaterialPropertyName>("L_BD")),
   _L_CD_name(getParam<MaterialPropertyName>("L_CD")),
   //Derivative with respect to xB
   _dL_BB_xB_name(getParam<MaterialPropertyName>("dL_BB_xB")),
   _dL_CC_xB_name(getParam<MaterialPropertyName>("dL_CC_xB")),
   _dL_DD_xB_name(getParam<MaterialPropertyName>("dL_DD_xB")),
   _dL_BC_xB_name(getParam<MaterialPropertyName>("dL_BC_xB")),
   _dL_BD_xB_name(getParam<MaterialPropertyName>("dL_BD_xB")),
   _dL_CD_xB_name(getParam<MaterialPropertyName>("dL_CD_xB")),
    //Derivative with respect to xC
   _dL_BB_xC_name(getParam<MaterialPropertyName>("dL_BB_xC")),
   _dL_CC_xC_name(getParam<MaterialPropertyName>("dL_CC_xC")),
   _dL_DD_xC_name(getParam<MaterialPropertyName>("dL_DD_xC")),
   _dL_BC_xC_name(getParam<MaterialPropertyName>("dL_BC_xC")),
   _dL_BD_xC_name(getParam<MaterialPropertyName>("dL_BD_xC")),
   _dL_CD_xC_name(getParam<MaterialPropertyName>("dL_CD_xC")),
   //Derivative with respect to xD
   _dL_BB_xD_name(getParam<MaterialPropertyName>("dL_BB_xD")),
   _dL_CC_xD_name(getParam<MaterialPropertyName>("dL_CC_xD")),
   _dL_DD_xD_name(getParam<MaterialPropertyName>("dL_DD_xD")),
   _dL_BC_xD_name(getParam<MaterialPropertyName>("dL_BC_xD")),
   _dL_BD_xD_name(getParam<MaterialPropertyName>("dL_BD_xD")),
   _dL_CD_xD_name(getParam<MaterialPropertyName>("dL_CD_xD")),
    //Declare properties
   _L_BB_val(declareProperty<Real>(_L_BB_name)),
   _L_CC_val(declareProperty<Real>(_L_CC_name)),
   _L_DD_val(declareProperty<Real>(_L_DD_name)),
   _L_BC_val(declareProperty<Real>(_L_BC_name)),
   _L_BD_val(declareProperty<Real>(_L_BD_name)),
   _L_CD_val(declareProperty<Real>(_L_CD_name)),
    //Derivatives wr.t xB
   _dL_BB_xB_val(declareProperty<Real>(_dL_BB_xB_name)),
   _dL_CC_xB_val(declareProperty<Real>(_dL_CC_xB_name)),
   _dL_DD_xB_val(declareProperty<Real>(_dL_DD_xB_name)),
   _dL_BC_xB_val(declareProperty<Real>(_dL_BC_xB_name)),
   _dL_BD_xB_val(declareProperty<Real>(_dL_BD_xB_name)),
   _dL_CD_xB_val(declareProperty<Real>(_dL_CD_xB_name)),
   //Derivatives wrt xC
   _dL_BB_xC_val(declareProperty<Real>(_dL_BB_xC_name)),
   _dL_CC_xC_val(declareProperty<Real>(_dL_CC_xC_name)),
   _dL_DD_xC_val(declareProperty<Real>(_dL_DD_xC_name)),
   _dL_BC_xC_val(declareProperty<Real>(_dL_BC_xC_name)),
   _dL_BD_xC_val(declareProperty<Real>(_dL_BD_xC_name)),
   _dL_CD_xC_val(declareProperty<Real>(_dL_CD_xC_name)),
   //Derivatives wrt xD
   _dL_BB_xD_val(declareProperty<Real>(_dL_BB_xD_name)),
   _dL_CC_xD_val(declareProperty<Real>(_dL_CC_xD_name)),
   _dL_DD_xD_val(declareProperty<Real>(_dL_DD_xD_name)),
   _dL_BC_xD_val(declareProperty<Real>(_dL_BC_xD_name)),
   _dL_BD_xD_val(declareProperty<Real>(_dL_BD_xD_name)),
   _dL_CD_xD_val(declareProperty<Real>(_dL_CD_xD_name)),
   //Each property is coupled to
   _xB(coupledValue("xB")),
   _xC(coupledValue("xC")),
   _xD(coupledValue("xD")),
   _table_object(getUserObject<QuaternaryMobilityData>("table_object"))
{
}

void
QuaternaryKineticMaterial::computeQpProperties()
{
    //Diagonal terms of the matrix
    _L_BB_val[_qp] = _table_object.L_BB(_xB[_qp], _xC[_qp], _xD[_qp]);
    _L_CC_val[_qp] = _table_object.L_CC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _L_DD_val[_qp] = _table_object.L_DD(_xB[_qp], _xC[_qp], _xD[_qp]);

    //off-diagonal terms
    _L_BC_val[_qp] = _table_object.L_BC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _L_BD_val[_qp] = _table_object.L_BD(_xB[_qp], _xC[_qp], _xD[_qp]);
    _L_CD_val[_qp] = _table_object.L_CD(_xB[_qp], _xC[_qp], _xD[_qp]);

    //Derivative with respect to xB
    _dL_BB_xB_val[_qp] = _table_object.dL_BB_xB(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CC_xB_val[_qp] = _table_object.dL_CC_xB(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_DD_xB_val[_qp] = _table_object.dL_DD_xB(_xB[_qp], _xC[_qp], _xD[_qp]);

    _dL_BC_xB_val[_qp] = _table_object.dL_BC_xB(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_BD_xB_val[_qp] = _table_object.dL_BD_xB(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CD_xB_val[_qp] = _table_object.dL_CD_xB(_xB[_qp], _xC[_qp], _xD[_qp]);

    //Derivative with respect to xC
    _dL_BB_xC_val[_qp] = _table_object.dL_BB_xC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CC_xC_val[_qp] = _table_object.dL_CC_xC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_DD_xC_val[_qp] = _table_object.dL_DD_xC(_xB[_qp], _xC[_qp], _xD[_qp]);

    _dL_BC_xC_val[_qp] = _table_object.dL_BC_xC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_BD_xC_val[_qp] = _table_object.dL_BD_xC(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CD_xC_val[_qp] = _table_object.dL_CD_xC(_xB[_qp], _xC[_qp], _xD[_qp]);

    //Derivative with respect to xD
    _dL_BB_xD_val[_qp] = _table_object.dL_BB_xD(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CC_xD_val[_qp] = _table_object.dL_CC_xD(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_DD_xD_val[_qp] = _table_object.dL_DD_xD(_xB[_qp], _xC[_qp], _xD[_qp]);

    _dL_BC_xD_val[_qp] = _table_object.dL_BC_xD(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_BD_xD_val[_qp] = _table_object.dL_BD_xD(_xB[_qp], _xC[_qp], _xD[_qp]);
    _dL_CD_xD_val[_qp] = _table_object.dL_CD_xD(_xB[_qp], _xC[_qp], _xD[_qp]);
}
