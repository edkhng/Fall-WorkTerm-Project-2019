#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4RunManager.hh"


///////////////////////////////////////////////////////////////////////////////

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* primary)
	: fPrimaryAction(primary)
{
  //Setup a command directory for detector controls with guidance
  fEventDir = new G4UIdirectory("/Ex/");
  fEventDir->SetGuidance("set eventid");

  // set event id
  fEventIDcmd = new G4UIcmdWithAnInteger("/Ex/eventID", this);
  fEventIDcmd->SetGuidance("Set Event ID");
  fEventIDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fEventIDcmd->SetToBeBroadcasted(false);

  // set prinary particle energy
  // fEnergycmd = new G4UIcmdWithADoubleAndUnit("/Ex/energy", this);
  // fEnergycmd->SetGuidance("Set Total Energy");
  // fEnergycmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  // fEnergycmd->SetToBeBroadcasted(false);

}

///////////////////////////////////////////////////////////////////////////////

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fEventIDcmd;
  // delete fEnergycmd;
  delete fEventDir;
}

///////////////////////////////////////////////////////////////////////////////

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if (command == fEventIDcmd) {
		fPrimaryAction->SetEventID(fEventIDcmd->GetNewIntValue(newValue));
	}
	// else if (command == fEnergycmd) {
	// 	fPrimaryAction->SetParticleEnergy(fEnergycmd->GetNewDoubleValue(newValue));
	// }
}
