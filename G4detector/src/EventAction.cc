#include "EventAction.hh"
#include "PMTSD.hh"
#include "PMTHit.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////

EventAction::EventAction()
 : G4UserEventAction(),
   fPMTHCID(-1)
{}

///////////////////////////////////////////////////////////////////////////////

EventAction::~EventAction()
{}

///////////////////////////////////////////////////////////////////////////////

PMTHitsCollection* EventAction::GetHitsCollection(G4int HCID,
	const G4Event* event) const
{
	auto hitsCollection
		= static_cast<PMTHitsCollection*>(
			event->GetHCofThisEvent()->GetHC(HCID));

	return hitsCollection;
}

///////////////////////////////////////////////////////////////////////////////

void EventAction::BeginOfEventAction(const G4Event* event)
{
	// Get hits collections IDs (only once)
	if (fPMTHCID == -1) {
		fPMTHCID
			= G4SDManager::GetSDMpointer()->GetCollectionID("PMTHitColl");
	}
}

///////////////////////////////////////////////////////////////////////////////

void EventAction::EndOfEventAction(const G4Event* event)
{
	// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();

	// Get hits collections
	auto hc = GetHitsCollection(fPMTHCID, event);
	if (!hc) return;
	auto nhit = hc->GetSize();
	G4int eid = event->GetEventID();
	analysisManager->FillNtupleIColumn(8, eid);

	G4cout << "End of Event, Event ID:" << eid << G4endl;
	G4cout<< "Number of Hit: " <<nhit<<G4endl;

	for (unsigned long i = 0; i < nhit; ++i) {
		// Fill Ntuple
		auto hit = static_cast<PMTHit*>(hc->GetHit(i));
		auto Time = hit->GetTime();
		auto Energy = hit->GetEnergy();
		auto WorldPos = hit->GetWorldPos();
		auto layerID = hit->GetLayerID();
		auto columnID = hit->GetColumnID();
		auto cellID = hit->GetCellID();
		
		analysisManager->FillNtupleIColumn(0, layerID);
		analysisManager->FillNtupleIColumn(1, columnID);
		analysisManager->FillNtupleIColumn(2, cellID);
		analysisManager->FillNtupleDColumn(3, Time);

		analysisManager->FillNtupleDColumn(4, WorldPos.x());
		analysisManager->FillNtupleDColumn(5, WorldPos.y());
		analysisManager->FillNtupleDColumn(6, WorldPos.z());
		analysisManager->FillNtupleDColumn(7, Energy);

		analysisManager->AddNtupleRow();
	}
}

///////////////////////////////////////////////////////////////////////////////
