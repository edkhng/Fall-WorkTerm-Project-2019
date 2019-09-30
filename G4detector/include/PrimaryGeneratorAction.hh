#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#include "PrimaryGeneratorMessenger.hh"

class PrimaryGenerator;
class G4Event;
class G4Box;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();    
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         

	void SetEventID(G4int ID) { eventID = ID; }
	void SetParticleEnergy(G4double Energy) { energy = Energy; }
	void SetParticleType(G4int Type) { type = Type; }
  
  private:
	PrimaryGenerator* fPrimaryGenerator;
    G4Box* worldBox;
	G4int eventID;
	G4int energy;
	G4int type;
	PrimaryGeneratorMessenger* fPrimaryGeneratorMessenger;
};


#endif
