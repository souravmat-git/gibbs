#include "SubstitutionalFreeEnergyMaterial.h"
#include "Material.h"
#include "cmath"

registerMooseObject("gibbsApp", SubstitutionalFreeEnergyMaterial);

//template <>
InputParameters
SubstitutionalFreeEnergyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<MaterialPropertyName>("free_energy_phase","Free energy of the phase");
  params.addRequiredParam<MaterialPropertyName>("B_diff_pot","First derivative w.r.t comp. B");
  params.addRequiredParam<MaterialPropertyName>("second_derivative", "Second derivative w.r.t B");
  params.addRequiredCoupledVar("mole_fraction_B","Mole fraction of B");
  params.addParam<Real>("char_energy", 1e9,"characteristic energy (J/m^3)");
  params.addParam<Real>("molar_volume", 1e-5, "Molar volume of the phase");
  params.addParam<Real>("gas_constant",8.314,"gas constant (J/molK)");
  params.addParam<Real>("temp",298,"temperature (K)");
  params.addRequiredParam<Real>("GHSER_A", "GHSER value of component A at 298K and 1atm");
  params.addRequiredParam<Real>("GHSER_B", "GHSER value of component B at 298K and 1atm");
  params.addRequiredParam<std::vector<std::string>>("RK_exp_terms", "Redlich-Kister expansion terms");
  params.addRequiredParam<std::vector<Real>>("RK_exp_values", "values of the expansion");
  params.addClassDescription("Computes the free energy and its derivatives of a substitutional solid solution");
  return params;
}

SubstitutionalFreeEnergyMaterial::SubstitutionalFreeEnergyMaterial(const InputParameters & parameters)
  : Material(parameters),
  _Ec(getParam<Real>("char_energy")),
  _Vm(getParam<Real>("molar_volume")),
  _R(getParam<Real>("gas_constant")),
  _T(getParam<Real>("temp")),
  _GHSERA(getParam<Real>("GHSER_A")),
  _GHSERB(getParam<Real>("GHSER_B")),
  _RK_names(getParam<std::vector<std::string>>("RK_exp_terms")),
  _RK_val(getParam<std::vector<Real>>("RK_exp_values")),
  _xB(coupledValue("mole_fraction_B")),
  _free_energy_name(getParam<MaterialPropertyName>("free_energy_phase")),
  _B_diff_pot_name(getParam<MaterialPropertyName>("B_diff_pot")),
  _chi_name(getParam<MaterialPropertyName>("second_derivative")),
  _free_energy(declareProperty<Real>(_free_energy_name)),
 _B_diff_pot(declareProperty<Real>(_B_diff_pot_name)),
  _chi(declareProperty<Real>(_chi_name))
{
  //See GenericFunctionMaterial.C for implementation

  if(_RK_names.size()!= _RK_val.size())

     mooseError("Redlich-Kister expansion terms and values must be equal");
}

void SubstitutionalFreeEnergyMaterial::computeQpProperties()
{
   _G_ref = (1.0 -_xB[_qp]) *_GHSERA + _xB[_qp]* _GHSERB;
   _G_conf = _R * _T * (_xB[_qp] * log(_xB[_qp]) + (1.0-_xB[_qp])* log(1.0 -_xB[_qp]));

    //Interaction energy
   _G_ex = (1.0 -_xB[_qp])*_xB[_qp]*(_RK_val[0] + _RK_val[1]*(1.0-2.0*_xB[_qp])
                                               + _RK_val[2]*std::pow((1.0-2.0*_xB[_qp]),2.0)
                                               + _RK_val[3]*std::pow((1.0-2.0*_xB[_qp]),3.0));

  _free_energy[_qp] = ((_G_ref + _G_conf + _G_ex)/_Vm)/_Ec ;

  //Compute the first derivative w.r.t component B.....

    _dG_ref = _GHSERB - _GHSERA;
    _dG_conf = _R* _T*(log(_xB[_qp]) - log(1.0 -_xB[_qp]));

   //Interaction energy can be split into two parts
   _dG_ex1 = (1.0-2.0*_xB[_qp])*(_RK_val[0]  + _RK_val[1]*(1.0-2.0*_xB[_qp])
                                             + _RK_val[2]*std::pow((1.0-2.0*_xB[_qp]),2.0)
                                             + _RK_val[3]*std::pow((1.0-2.0*_xB[_qp]),3.0));

   _dG_ex2 = -(1.0-_xB[_qp])*_xB[_qp]*2.0*(_RK_val[1]
                                           + 2.0*_RK_val[2] * (1.0-2.0*_xB[_qp])
                                           + 3.0*_RK_val[3] * std::pow((1.0-2.0*_xB[_qp]),2.0));

   _B_diff_pot[_qp] = ((_dG_ref + _dG_conf + _dG_ex1 + _dG_ex2)/_Vm)/_Ec ;

  //Compute the second derivatives....

  _d2G_ref = 0.0;

  _d2G_conf = _R* _T*(1.0/(_xB[_qp] * (1.0 -_xB[_qp])));

  _d2G_ex1 = -2.0*(_RK_val[0] + _RK_val[1]*(1.0-2.0*_xB[_qp])
                              + _RK_val[2]*std::pow((1.0-2.0*_xB[_qp]),2.0)
                              + _RK_val[3]*std::pow((1.0-2.0*_xB[_qp]),3.0))
             -2.0*(1.0 - 2.0*_xB[_qp])*(_RK_val[1] + 2.0*_RK_val[2]*(1.0 - 2.0*_xB[_qp])
                   + 3.0*_RK_val[3]*std::pow((1.0-2.0*_xB[_qp]),2.0));

  //dG_ex2 is calculated in three parts and finally summed
  _d2G_ex2 = -2.0*(1.0 -2.0*_xB[_qp])*(_RK_val[1]
                                       + 2.0* _RK_val[2]* (1.0-2.0*_xB[_qp])
                                       + 3.0* _RK_val[3]*std::pow((1.0 -2.0*_xB[_qp]),2.0));

  _d2G_ex3 = 4.0*(1.0-_xB[_qp])*_xB[_qp]*(2.0*_RK_val[2] + 6.0*_RK_val[3]*(1.0-2.0*_xB[_qp]));

  _chi[_qp] = ((_d2G_ref + _d2G_conf + _d2G_ex1 + _d2G_ex2 + _d2G_ex3)/_Vm)/_Ec;

}
