#include "PrimaryGenerator.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include<math.h>

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

		// CLHEP::HepRandomEngine* theEngine = CLHEP::HepRandom::getTheEngine(); // Random engin
		// CLHEP::RandExponential* Randexp = new CLHEP::RandExponential(theEngine);

		for (int i = 0; i < 5000000; i++) {
			// G4double phi = G4UniformRand() * 360 * deg; //uniformly distributed
			// G4double theta = Randexp->shoot() / 3.1415926 * 0.5 * deg; //exponential distribution

			G4double phi = G4UniformRand() * 2 * std::acos(-1);
			G4double theta = std::acos(1 - 0.05 * G4UniformRand());

			G4double uz = - std::cos(theta);
			G4double ux = std::sin(theta) * std::cos(phi);
			G4double uy = std::sin(theta) * std::sin(phi);
			G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);

      G4double x = ux*40/sqrt(2969) + uz*37/sqrt(2969);
      G4double z = -ux*37/sqrt(2969) + uz*40/sqrt(2969);

      particle1->SetMomentumDirection(G4ThreeVector(x, uy, z));
			particle1->SetPolarization(G4ThreeVector(ux, uy, uz));
			//particle1->SetMomentumDirection(G4ThreeVector(0, 0, -1));
			//particle1->SetPolarization(G4ThreeVector(0, 0, -1));
			particle1->SetKineticEnergy(3.061*eV);
			vertexA->SetPrimary(particle1);
		}
		event->AddPrimaryVertex(vertexA);
}
