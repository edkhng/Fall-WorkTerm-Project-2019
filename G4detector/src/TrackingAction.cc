#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	if (aTrack->GetParentID() == 0)
	{	// true: only stroe trajectory information for the primary particle
		fpTrackingManager->SetStoreTrajectory(false);
	}
	else
	{
		fpTrackingManager->SetStoreTrajectory(false);
		// used to monitor the process
		//if (aTrack->GetParentID() == 1 && aTrack->GetCreatorProcess()->GetProcessName() == "Decay") {
		//	G4ThreeVector pos = aTrack->GetVertexPosition();
		//	G4cout << "delta_x: " << pos.x() + 10000 << "  ";
		//	G4cout << "delta_y: " << pos.y() + 10000 << "  ";
		//	G4cout << "delta_z: " << pos.z() + 125000 << "  ";
		//	G4cout << aTrack->GetParticleDefinition()->GetParticleName() << G4endl;
		//}
	}
}


