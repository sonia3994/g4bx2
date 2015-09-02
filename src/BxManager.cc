//Created by D. Franco
//Revised by A. Caminata and S. Marcocci, Sept. 2014
#include "G4RunManager.hh"
#include "BxLogger.hh"
#include "BxManagerMessenger.hh"
#include "Randomize.hh"
#include <time.h>
#include <iostream>
#include "BxManager.hh"
#include "BxIO.hh"

//---------------------------------------------------------------------------//
BxManager* BxManager::Manager = 0;

BxManager::BxManager() {
  Manager         = this;   
  fBxMessenger    = new BxManagerMessenger(this);
}

BxManager* BxManager::Get() {
  if (!Manager)
    Manager = new BxManager();
  return Manager;
}

BxManager::~BxManager()
{
  delete fBxMessenger;
}
