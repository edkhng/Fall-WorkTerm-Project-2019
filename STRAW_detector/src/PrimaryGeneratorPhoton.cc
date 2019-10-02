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
	// Uniform light source
	G4ParticleDefinition* OpPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");
	G4double phi = G4UniformRand() * 360 * deg;
	G4double theta = G4UniformRand() * 180 * deg;
	G4double uz = std::cos(theta);
	G4double ux = std::sin(theta) * std::cos(phi);
	G4double uy = std::sin(theta) * std::sin(phi);
	G4PrimaryParticle* particle = new G4primaryParticle(OpPhoton);
		vertexB->AddPrimary(particle);
	event->AddprimaryVertex(vertexB)

}
