#include "gibbsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
gibbsApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  //params.set<bool>("use_legacy_material_output") = false;

  return params;
}

gibbsApp::gibbsApp(InputParameters parameters) : MooseApp(parameters)
{
  gibbsApp::registerAll(_factory, _action_factory, _syntax);
}

gibbsApp::~gibbsApp() {}

void
gibbsApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"gibbsApp"});
  Registry::registerActionsTo(af, {"gibbsApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
gibbsApp::registerApps()
{
  registerApp(gibbsApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
gibbsApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  gibbsApp::registerAll(f, af, s);
}
extern "C" void
gibbsApp__registerApps()
{
  gibbsApp::registerApps();
}
