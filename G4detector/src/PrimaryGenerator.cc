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
	if (ftype == 1) {
		//vertex A
		G4double timeA = 0 * s;
		//
		G4PrimaryVertex* vertexA = new G4PrimaryVertex(fpos, timeA);

		//Charged lepton at vertex A
		//
		G4ParticleDefinition* particleDefinition
			= G4ParticleTable::GetParticleTable()->FindParticle("tau-");
		G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
		particle1->SetMomentumDirection(G4ThreeVector(0, 0, 1));
		//particle1->SetMomentumDirection(G4ThreeVector(0, 1, 0));
		particle1->SetKineticEnergy(fenergy);
		//
		vertexA->SetPrimary(particle1);
		//G4cout << "tau: " << fenergy / TeV << G4endl;
		event->AddPrimaryVertex(vertexA);
	}
	if (ftype == 2) {
		//vertex A
		G4double timeA = 0 * s;
		//
		G4PrimaryVertex* vertexA = new G4PrimaryVertex(fpos, timeA);

		//Charged lepton at vertex A
		//
		G4ParticleDefinition* particleDefinition
			= G4ParticleTable::GetParticleTable()->FindParticle("e-");
		G4PrimaryParticle* particle1 = new G4PrimaryParticle(particleDefinition);
		particle1->SetMomentumDirection(G4ThreeVector(0, 0, 1));
		//particle1->SetMomentumDirection(G4ThreeVector(0, 1, 0));
		particle1->SetKineticEnergy(fenergy);
		//
		vertexA->SetPrimary(particle1);
		//G4cout << "e: " << fenergy / TeV << G4endl;
		event->AddPrimaryVertex(vertexA);
	}


	else if (ftype == 0) {
		//vertex B
		G4double timeB = 0 * s;
		G4PrimaryVertex* vertexB = new G4PrimaryVertex(fpos, timeB);

		G4ParticleDefinition* Pi_nu = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
		G4ParticleDefinition* Pi_pu = G4ParticleTable::GetParticleTable()->FindParticle("pi+");
		G4ParticleDefinition* Pi_ne = G4ParticleTable::GetParticleTable()->FindParticle("pi-");

		CLHEP::HepRandomEngine* theEngine = CLHEP::HepRandom::getTheEngine(); // Random engin
		CLHEP::RandExponential* Randexp = new CLHEP::RandExponential(theEngine);
		G4cout << "hadronic:" << fenergy / TeV << G4endl;

		//pions at vertex B
		for (int i = 0; i < 20; i++) {
			G4double phi = G4UniformRand() * 360 * deg; //uniformly distributed
			G4double theta = Randexp->shoot() / 3.1415926 * 0.36 * deg; //exponential distribution

			// calculate the momentum direction
			G4double uz = std::cos(theta);
			G4double ux = std::sin(theta) * std::cos(phi);
			G4double uy = std::sin(theta) * std::sin(phi);
			//G4double uy = std::cos(theta);
			//G4double ux = std::sin(theta) * std::cos(phi);
			//G4double uz = std::sin(theta) * std::sin(phi);
			G4PrimaryParticle* particle2;

			G4double flag = G4UniformRand(); // Assign random charge to pions
			if (flag < 0.333333) {
				particle2 = new G4PrimaryParticle(Pi_nu);
				G4cout << "Pi0:" << G4endl;
			}
			else if (flag < 0.666667) {
				particle2 = new G4PrimaryParticle(Pi_pu);
				G4cout << "Pi+:" << G4endl;
			}
			else {
				particle2 = new G4PrimaryParticle(Pi_ne);
				G4cout << "Pi-:" << G4endl;
			}

			//G4cout << "Phi: " << phi << "   Theta: " << theta << G4endl;
			particle2->SetMomentumDirection(G4ThreeVector(ux, uy, uz));
			particle2->SetKineticEnergy(fenergy/20);
			vertexB->SetPrimary(particle2);

		}
		event->AddPrimaryVertex(vertexB);

	}

	/* Uniform light source
	G4ParticleDefinition* OpPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
	G4double phi = G4UniformRand() * 360 * deg;
	G4double theta = G4UniformRand() * 180 * deg;
	G4double uz = std::cos(theta);
	G4double ux = std::sin(theta) * std::cos(phi);
	G4double uy = std::sin(theta) * std::sin(phi);
	G4PrimaryParticle* particle = new G4primaryParticle(OpPhoton);
		vertexB->AddPrimary(particle);
	event->AddprimaryVertex(vertexB)
	*/

}
