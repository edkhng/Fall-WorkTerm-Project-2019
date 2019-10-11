#include "PMTHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


#include <iomanip>

G4ThreadLocal G4Allocator<PMTHit>* PMTHitAllocator = 0;

///////////////////////////////////////////////////////////////////////////////

PMTHit::PMTHit()
  : G4VHit(),
    fTime(0.)
{}

PMTHit::PMTHit(const PMTHit &right)
: G4VHit(),
  fTime(right.fTime)
{}

///////////////////////////////////////////////////////////////////////////////

PMTHit::~PMTHit()
{}

///////////////////////////////////////////////////////////////////////////////

const PMTHit& PMTHit::operator=(const PMTHit &right)
{
  fTime = right.fTime;
  return *this;
}

///////////////////////////////////////////////////////////////////////////////

G4bool PMTHit::operator==(const PMTHit &right) const
{
  return ( this == &right ) ? true : false;
}

///////////////////////////////////////////////////////////////////////////////

// void PMTHit::Draw()
// {
// 	// /vis/scene/add/hits
// 	auto visManager = G4VVisManager::GetConcreteInstance();
// 	if (!visManager) return;
//
// 	G4Circle circle(fWorldPos);
// 	circle.SetScreenSize(2);
// 	circle.SetFillStyle(G4Circle::filled);
// 	G4Colour colour(1., 0., 1.);
// 	G4VisAttributes attribs(colour);
// 	circle.SetVisAttributes(attribs);
// 	visManager->Draw(circle);
// }

///////////////////////////////////////////////////////////////////////////////

const std::map<G4String, G4AttDef>* PMTHit::GetAttDefs() const
{
	G4bool isNew;
	auto store = G4AttDefStore::GetInstance("PMTHit", isNew);

	if (isNew) {
		(*store)["HitType"]
			= G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");

		(*store)["Time"]
			= G4AttDef("Time", "Time", "Physics", "G4BestUnit", "G4double");
	}

	return store;
}

///////////////////////////////////////////////////////////////////////////////

std::vector<G4AttValue>* PMTHit::CreateAttValues() const
{
	auto values = new std::vector<G4AttValue>;

	values
		->push_back(G4AttValue("HitType", "DriftChamberHit", ""));
	values
		->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));

	return values;
}

///////////////////////////////////////////////////////////////////////////////

// void PMTHit::Print()
// {
// 	G4cout << "  Layer["
// 		<< fLayerID << "]["
// 		<< fColumnID << "]["
// 		<< fCellID << "]"
// 		<< " : time " << fTime / ns << G4endl;
// }

///////////////////////////////////////////////////////////////////////////////
