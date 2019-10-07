#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGenerator.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


///////////////////////////////////////////////////////////////////////////////

PrimaryGeneratorAction::PrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fPrimaryGenerator(0),
	eventID(0),
	energy(0),
	type(0)
{
  fPrimaryGenerator = new PrimaryGenerator();
  fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);
}

///////////////////////////////////////////////////////////////////////////////

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fPrimaryGenerator;
}

///////////////////////////////////////////////////////////////////////////////

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  anEvent->SetEventID(eventID);
  eventID++;
  CLHEP::HepRandom::setTheSeed(eventID);

  // store world volume information
  G4double worldHalfSizeX = 0;
  G4double worldHalfSizeY = 0;
  G4double worldHalfSizeZ = 0;

  // get world volume from G4LogicalVolumeStore
  G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("WorldLV");
  
  if ( worldLV ) {
	  worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
	worldHalfSizeX = worldBox->GetXHalfLength();
	worldHalfSizeY = worldBox->GetYHalfLength();
    worldHalfSizeZ = worldBox->GetZHalfLength();
  }  

  // Set gun position, energy, type
  fPrimaryGenerator->SetParticlePosition(G4ThreeVector(	18.5 * m, 0 * m, 60 * m));
  fPrimaryGenerator->SetParticleEnergy(energy);
  fPrimaryGenerator->SetParticleType(type);
  fPrimaryGenerator->GeneratePrimaryVertex(anEvent);
}

///////////////////////////////////////////////////////////////////////////////