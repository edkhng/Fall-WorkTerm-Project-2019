#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

///////////////////////////////////////////////////////////////////////////////

PrimaryGenerator::PrimaryGenerator()
: G4VPrimaryGenerator(),
  fpos(0,0,0),
  fenergy(0),
  ftype(0)
{ }

///////////////////////////////////////////////////////////////////////////////

PrimaryGenerator::~PrimaryGenerator()
{ }

///////////////////////////////////////////////////////////////////////////////

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
		//vertex A
		G4double timeA = 0 * s;
		// 
		G4PrimaryVertex* vertexA = new G4PrimaryVertex(fpos, timeA);

		//Photon at A
		//
		G4ParticleDefinition* particleDefinition
			= G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
		for (int i = 0; i < 1; i++) {
			G4double phi = G4UniformRand() * 2 * 3.14159;
			G4double theta = std::acos(1 - 2 * G4UniformRand());

			G4double uz = std::cos(theta);
			G4double ux = std::sin(theta) * std::cos(phi);
			G4double uy = std::sin(theta) * std::sin(phi);
			G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);

			particle1->SetMomentumDirection(G4ThreeVector(ux, uy, uz));
			//particle1->SetPolarization(G4ThreeVector(ux, uy, uz));
			//particle1->SetMomentumDirection(G4ThreeVector(0, 0, -1));
			particle1->SetKineticEnergy(fenergy);
			vertexA->SetPrimary(particle1);
		}
		event->AddPrimaryVertex(vertexA);
}

