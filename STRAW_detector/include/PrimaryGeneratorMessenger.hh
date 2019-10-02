#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcommand;


/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

class PrimaryGeneratorMessenger : public G4UImessenger
{
  public:

    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
    virtual ~PrimaryGeneratorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);
  
  private:
	PrimaryGeneratorAction* fPrimaryAction;
	G4UIdirectory* fEventDir;
	G4UIcmdWithAnInteger* fEventIDcmd;
	G4UIcmdWithADoubleAndUnit* fEnergycmd;
	G4UIcmdWithAnInteger* fTypecmd;
};


#endif
