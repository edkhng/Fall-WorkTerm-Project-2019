#include "TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	if (aTrack->GetParentID() == 0)
	{	// true: only stroe trajectory information for the primary particle
		fpTrackingManager->SetStoreTrajectory(true);
	}
	else
	{
		fpTrackingManager->SetStoreTrajectory(false);
		// used to monitor the process
		if (aTrack->GetParentID() == 1 && aTrack->GetCreatorProcess()->GetProcessName() == "Decay" && aTrack->GetParticleDefinition()->GetParticleName() == "nu_tau")
		// if (aTrack->GetParentID() == 1 && aTrack->GetCreatorProcess()->GetProcessName() == "Decay" )
		{
			G4ThreeVector pos = aTrack->GetVertexPosition();
			G4cout << "delta_x: " << pos.x() + 10000 << " [mm] ;";
			G4cout << "delta_y: " << pos.y() + 20000 << " [mm] ";;
			G4cout << "delta_z: " << pos.z() + 62500 << " [mm] ;";
			G4cout << "decay_time: " << (pos.z() + 62500)/299.792458 << " [ns] ";
			G4cout << aTrack->GetParticleDefinition()->GetParticleName() << G4endl;
		}
	}
}
