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
  fLayerID(-1), fColumnID(-1), fCellID(-1), fTime(0.), fEnergy(0), fWorldPos(0)
{}

PMTHit::PMTHit(G4int layerID, G4int columnID, G4int cellID)
: G4VHit(),
  fLayerID(layerID), fColumnID(columnID), fCellID(cellID),
  fTime(0.), fEnergy(0),
  fWorldPos(0)
{}

PMTHit::PMTHit(const PMTHit &right)
: G4VHit(),
  fLayerID(right.fLayerID),
  fColumnID(right.fColumnID),
  fCellID(right.fCellID),
  fTime(right.fTime),
 fEnergy(right.fEnergy),
  fWorldPos(right.fWorldPos)
{}

///////////////////////////////////////////////////////////////////////////////

PMTHit::~PMTHit()
{}

///////////////////////////////////////////////////////////////////////////////

const PMTHit& PMTHit::operator=(const PMTHit &right)
{
  fLayerID = right.fLayerID;
  fColumnID = right.fColumnID;
  fCellID = right.fCellID;
  fTime = right.fTime;
  fEnergy = right.fEnergy;
  fWorldPos = right.fWorldPos;
  return *this;
}

///////////////////////////////////////////////////////////////////////////////

G4bool PMTHit::operator==(const PMTHit &right) const
{
  return ( this == &right ) ? true : false;
}

///////////////////////////////////////////////////////////////////////////////

void PMTHit::Draw()
{
	// /vis/scene/add/hits
	auto visManager = G4VVisManager::GetConcreteInstance();
	if (!visManager) return;

	G4Circle circle(fWorldPos);
	circle.SetScreenSize(2);
	circle.SetFillStyle(G4Circle::filled);
	G4Colour colour(1., 0., 1.);
	G4VisAttributes attribs(colour);
	circle.SetVisAttributes(attribs);
	visManager->Draw(circle);
}

///////////////////////////////////////////////////////////////////////////////

const std::map<G4String, G4AttDef>* PMTHit::GetAttDefs() const
{
	G4bool isNew;
	auto store = G4AttDefStore::GetInstance("PMTHit", isNew);

	if (isNew) {
		(*store)["HitType"]
			= G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");

		(*store)["Layer"]
			= G4AttDef("Layer", "Layer ID", "Physics", "", "G4int");

		(*store)["Column"]
			= G4AttDef("Column", "Column ID", "Physics", "", "G4int");

		(*store)["Cell"]
			= G4AttDef("Cell", "Cell ID", "Physics", "", "G4int");

		(*store)["Time"]
			= G4AttDef("Time", "Time", "Physics", "G4BestUnit", "G4double");

		(*store)["Pos"]
			= G4AttDef("Pos", "Position", "Physics", "G4BestUnit", "G4ThreeVector");
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
		->push_back(G4AttValue("Layer", G4UIcommand::ConvertToString(fLayerID), ""));
	values
		->push_back(G4AttValue("Column", G4UIcommand::ConvertToString(fLayerID), ""));
	values
		->push_back(G4AttValue("Cell", G4UIcommand::ConvertToString(fLayerID), ""));
	values
		->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));
	values
		->push_back(G4AttValue("Pos", G4BestUnit(fWorldPos, "Length"), ""));

	return values;
}

///////////////////////////////////////////////////////////////////////////////

void PMTHit::Print()
{
	G4cout << "  Layer["
		<< fLayerID << "]["
		<< fColumnID << "]["
		<< fCellID << "]"
		<< " : time " << fTime / ns << G4endl;
}

///////////////////////////////////////////////////////////////////////////////