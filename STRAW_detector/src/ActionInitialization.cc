#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"
#include "TrackingAction.hh"


///////////////////////////////////////////////////////////////////////////////

ActionInitialization::ActionInitialization()
 : G4VUserActionInitialization()
{}

///////////////////////////////////////////////////////////////////////////////

ActionInitialization::~ActionInitialization()
{}

///////////////////////////////////////////////////////////////////////////////

void ActionInitialization::BuildForMaster() const
{

}

///////////////////////////////////////////////////////////////////////////////

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction);
  SetUserAction(new RunAction());
  EventAction* eventAction = new EventAction();

  SetUserAction(eventAction); 
  SetUserAction(new StackingAction());
  SetUserAction(new TrackingAction());
}  

