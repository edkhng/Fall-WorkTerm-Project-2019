#include "StackingAction.hh"

#include "G4VProcess.hh"

#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

///////////////////////////////////////////////////////////////////////////////

StackingAction::StackingAction()
  : G4UserStackingAction(),
    fCerenkovCounter(0)
{}

///////////////////////////////////////////////////////////////////////////////

StackingAction::~StackingAction()
{}

///////////////////////////////////////////////////////////////////////////////

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
		// particle is cerenkov photon
        fCerenkovCounter++;
    }
  }
  return fUrgent;
}

///////////////////////////////////////////////////////////////////////////////

void StackingAction::NewStage()
{
  G4cout << "Number of Cerenkov photons produced in this event : "
         << fCerenkovCounter << G4endl;
}

///////////////////////////////////////////////////////////////////////////////

void StackingAction::PrepareNewEvent()
{
  fCerenkovCounter = 0;
}

///////////////////////////////////////////////////////////////////////////////
