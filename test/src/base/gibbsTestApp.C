//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "gibbsTestApp.h"
#include "gibbsApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
gibbsTestApp::validParams()
{
  InputParameters params = gibbsApp::validParams();
  return params;
}

gibbsTestApp::gibbsTestApp(InputParameters parameters) : MooseApp(parameters)
{
  gibbsTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

gibbsTestApp::~gibbsTestApp() {}

void
gibbsTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  gibbsApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"gibbsTestApp"});
    Registry::registerActionsTo(af, {"gibbsTestApp"});
  }
}

void
gibbsTestApp::registerApps()
{
  registerApp(gibbsApp);
  registerApp(gibbsTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
gibbsTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  gibbsTestApp::registerAll(f, af, s);
}
extern "C" void
gibbsTestApp__registerApps()
{
  gibbsTestApp::registerApps();
}
