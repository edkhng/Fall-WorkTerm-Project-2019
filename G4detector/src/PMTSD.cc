#include "PMTSD.hh"
#include "PMTHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

///////////////////////////////////////////////////////////////////////////////

PMTSD::PMTSD(G4String name)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr),
   fHCID(-1)
{
  collectionName.insert("PMTHitColl");
}

///////////////////////////////////////////////////////////////////////////////

PMTSD::~PMTSD()
{}

///////////////////////////////////////////////////////////////////////////////

void PMTSD::Initialize(G4HCofThisEvent* HCE)
{
  // Create hits collection
  fHitsCollection 
    = new PMTHitsCollection(SensitiveDetectorName, collectionName[0]);

  if (fHCID < 0) {
	  fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  
  // Add this collection in hce
  HCE->AddHitsCollection( fHCID, fHitsCollection ); 
}

///////////////////////////////////////////////////////////////////////////////

G4bool PMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{  
	auto preStepPoint = step->GetPreStepPoint();
	auto track = step->GetTrack();

	// get cell, column, and layer number
	auto touchable = step->GetPreStepPoint()->GetTouchable();
	auto cellNumber = touchable->GetReplicaNumber(1);
	auto columnNumber = touchable->GetReplicaNumber(2);
	auto layerNumber = touchable->GetReplicaNumber(3);

	auto hit = new PMTHit(layerNumber, columnNumber, cellNumber);
	hit->SetWorldPos(preStepPoint->GetPosition());
	hit->SetTime(preStepPoint->GetGlobalTime());
	hit->SetEnergy(preStepPoint->GetTotalEnergy());
	
	fHitsCollection->insert(hit);
	/*
	G4cout
		<< G4endl
		<< "Copy number" << layerNumber << G4endl
		<< "local position" << localPos << G4endl
		<< "time" << preStepPoint->GetGlobalTime() << G4endl;
	*/
	// kill the particle
	track->SetTrackStatus(fStopAndKill);
	return true;
}

///////////////////////////////////////////////////////////////////////////////

void PMTSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( print ) { 
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl 
       << "-------->Hits Collection: in this event they are " << nofHits 
       << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

///////////////////////////////////////////////////////////////////////////////
