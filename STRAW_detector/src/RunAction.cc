#include "RunAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
///////////////////////////////////////////////////////////////////////////////

RunAction::RunAction()
 : G4UserRunAction()
{
  G4RunManager::GetRunManager()->SetPrintProgress(1);
  // Create analysis manager
  // Analysis technology specified in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Creating ntuple
  G4cout<<"creat ntuple"<<G4endl;
  analysisManager->CreateNtuple("Ntuple", "hits");
  analysisManager->CreateNtupleDColumn("Time");       // column Id = 0

  analysisManager->FinishNtuple();

  //analysisManager->CreateNtupleDColumn(1, "hit_per_event");  // column Id = 0
  //analysisManager->FinishNtuple(1);
}

///////////////////////////////////////////////////////////////////////////////

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();
}

///////////////////////////////////////////////////////////////////////////////

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //auto anEvent = G4RunManager::GetRunManager()->GetCurrentEvent()->SetEventID(2);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile();
}

///////////////////////////////////////////////////////////////////////////////

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  //
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}

///////////////////////////////////////////////////////////////////////////////
