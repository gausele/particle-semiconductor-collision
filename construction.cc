#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction(){}
MyDetectorConstruction::~MyDetectorConstruction(){}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    G4NistManager *nist = G4NistManager::Instance();

    G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box *solidWorld = new G4Box("solidWorld", 0.5*m, 0.5*m, 0.5*m);
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld,
                                                     "physWorld", 0, false, 0, true);

    G4Element* Ga = nist->FindOrBuildElement("Ga");
    G4Element* N = nist->FindOrBuildElement("N");

    G4double density = 6.15*g/cm3;
    G4Material* GaN = new G4Material("GaN", density, 2);
    GaN->AddElement(Ga, 1);
    GaN->AddElement(N, 1);

    G4Box *solidTarget = new G4Box("solidTarget", 1.*cm, 1.*cm, 0.05*cm);
    G4LogicalVolume *logicTarget = new G4LogicalVolume(solidTarget, GaN, "logicTarget");    
    new G4PVPlacement(0, G4ThreeVector(0,0,5.*cm), logicTarget,
                      "physTarget", logicWorld, false, 0, true);


    G4VisAttributes *visTarget = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
    logicTarget->SetVisAttributes(visTarget);


    return physWorld;    
}

