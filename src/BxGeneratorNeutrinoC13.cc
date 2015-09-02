// --------------------------------------------------------------------------//
/** 
 * AUTHOR: I.Machulin
 * CONTACT: machulin@lngs.infn.it
 * Revised by A. Caminata and S. Marcocci, Sept. 2014
*/
// --------------------------------------------------------------------------//

#include "BxEventAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Neutron.hh"
#include "BxLogger.hh"
#include "BxDataCollection.hh"
#include "BxOutputVertex.hh"
#include "BxGeneratorNeutrinoC13.hh"


#include "CLHEP/Random/RandExponential.h"

namespace CLHEP {} 
using namespace CLHEP;

//---------------------------------------------------------------------------//


BxGeneratorNeutrinoC13::BxGeneratorNeutrinoC13(): BxVGenerator("BxGeneratorNeutrinoC13") {

  isFirstTime = true ;
  fNeutrinoType = -1;
  


  fVolumeFlag = false ;
  G4ThreeVector zero(0., 0., 0.) ;
  fSPSAng = new G4SPSAngDistribution ;
  fSPSPos = new G4SPSPosDistribution ;
  G4SPSRandomGenerator *RndGen = new G4SPSRandomGenerator;
  fSPSPos->SetBiasRndm(RndGen);
  fSPSAng->SetBiasRndm(RndGen);
  SetRadius0(0.);
  SetCentreCoords(zero);    
  fPosition = G4ThreeVector(0.,0.,0.) ;
  fParticleTable = G4ParticleTable::GetParticleTable();    

  fTheMessenger = new BxGeneratorNeutrinoC13Messenger(this);
 
}
//---------------------------------------------------------------------------//


BxGeneratorNeutrinoC13::~BxGeneratorNeutrinoC13()
{
  delete fTheMessenger;
  delete fSPSPos;
  delete fSPSAng;
  
}

//---------------------------------------------------------------------------//

void BxGeneratorNeutrinoC13::BxGeneratePrimaries(G4Event *event) {

  if(isFirstTime) {
    if(fNeutrinoType < 0) {
      BxLog(error) << "Set the neutrino + C13 source type"<<endlog;
      BxLog(fatal) << endlog;
    }  

    BxLog(routine) << "NeutrinoC13 source " << fNeutrinoType << endlog ;

    isFirstTime = false ;
  }
//  fNeutrinoType {FullN=0,FirstN=1,SecondN=2,FlatN=3};
//  Start of the  first event - nu + C13 - > e- + N13  

   if ((event->GetEventID())%2 == 0){

  if(fVolumeFlag)          
  {
  fSPSAng->SetAngDistType("iso");   
  fSPSAng->SetPosDistribution(fSPSPos); 
  fPosition = fSPSPos->GenerateOne();
  }  

             
    fParticle = fParticleTable->FindParticle(11);     // 11 for the electron
    
    G4double KinE=0;
    
    if ((fNeutrinoType == 0)||(fNeutrinoType == 1)) KinE = ShootEnergyElectronNuC13()/MeV;

    if (fNeutrinoType == 2) KinE = 0./MeV; 
    
    if (fNeutrinoType == 3) KinE = ShootEnergyElectronFlat()/MeV;

         
    G4double energy = KinE + fParticle->GetPDGMass();	
    
    fDirection = fSPSAng->GenerateOne();
    
    G4double pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();

    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
    if ((fNeutrinoType == 0)||(fNeutrinoType == 1)||(fNeutrinoType == 3))    
    {
    BxOutputVertex::Get()->SetIsotopeCoinc(0);
    }    
    
    
    fPositionFirst = fPosition ;   

    BxOutputVertex::Get()->SetEnergy(KinE);		
    BxOutputVertex::Get()->SetPosition(fPosition);		
    BxOutputVertex::Get()->SetDirection(fDirection);
    BxOutputVertex::Get()->SetTime(0);
    BxOutputVertex::Get()->SetPDG(fParticle->GetPDGEncoding());

   }


   //End of the first event generation

   else {
	   if(fNeutrinoType!=1){
		   // Second event - beta+ decay of N13

		   fPosition = fPositionFirst; 

		   G4double fDecayTime= 14.37 * 60.0 * 1000000000 ; //14.37 min in nsec: Decay time of N13
		   G4double delay = RandExponential::shoot(fDecayTime);


		   //    G4double delay = 200000 ; // The delay fixed to 200 mcs

		   //	theDataCollection->SetNeutronTime(delay);
		   //	theDataCollection->SetNeutronCaptureFlag(1);


		   fParticle = fParticleTable->FindParticle(-11); //positron     

		   G4double KinE = ShootEnergyPositronN13()/MeV;


		   G4double energy = KinE + fParticle->GetPDGMass();

		   fDirection = fSPSAng->GenerateOne();

		   G4double pmom = std::sqrt(pow(energy,2.) - pow(fParticle->GetPDGMass(),2.));
		   G4double px = pmom*fDirection.x();
		   G4double py = pmom*fDirection.y();
		   G4double pz = pmom*fDirection.z();


		   G4PrimaryVertex* vertex   = new G4PrimaryVertex(fPosition,0);
		   G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);

		   vertex->SetPrimary( particle );
		   event->AddPrimaryVertex( vertex );    



		   if ((fNeutrinoType == 0)||(fNeutrinoType == 3))    

		   {
			   BxOutputVertex::Get()->SetIsotopeCoinc(1);
		   } 

		   BxOutputVertex::Get()->SetEnergy(KinE);		
		   BxOutputVertex::Get()->SetPosition(fPosition);		
		   BxOutputVertex::Get()->SetDirection(fDirection);
		   BxOutputVertex::Get()->SetPDG(fParticle->GetPDGEncoding());
		   BxOutputVertex::Get()->SetNSequence(0);
		   BxOutputVertex::Get()->SetTime(delay);

	   }
   }

}

//-------------------------------------------------------------------------
//    Energy Shooter for electrons from nu+C13 -> e- + N13
//-------------------------------------------------------------------------

G4double BxGeneratorNeutrinoC13::ShootEnergyElectronNuC13() {

   G4double sum = 0.;
   G4double norma = 0.;
   G4double Probability [401];
   G4double EnergyBin [401];
   G4double E0 = 14.0; //end Total Energy of the spectre
   G4double scale = E0/400;
   G4double xE;
   G4double deltaX, x, y;

//calculated by A. Ianni
   G4double ElectronEnergySp [401] =    
 {   
0.00015145492627697300,
0.00029556796725553700,
0.00041440237019949900,
0.00053208997268673100,
0.00065416147162235000,
0.00078285659095058100,
0.00091947348005629000,
0.00106476340798130000,
0.00121930751287220000,
0.00138364670432791000,
0.00155808501172612000,
0.00174301610677081000,
0.00193879015908930000,
0.00214558496352623000,
0.00236378269884324000,
0.00259359830914266000,
0.00283515481829909000,
0.00308889957237249000,
0.00335489041741843000,
0.00363323978990324000,
0.00392445258768807000,
0.00422844322374562000,
0.00454538030623088000,
0.00487564502934270000,
0.00521914324772110000,
0.00557612531625113000,
0.00594682354317072000,
0.00633114144455952000,
0.00672942029497035000,
0.00714174043247372000,
0.00756799348277843000,
0.00800859163365911000,
0.00846342008753985000,
0.00893244507006154000,
0.00941617401258578000,
0.00991430044560414000,
0.01042689742307740000,
0.01095430287343770000,
0.01149616182648810000,
0.01205276759110070000,
0.01262421812110250000,
0.01321009388312610000,
0.01381103210432620000,
0.01442672605367770000,
0.01505676343960460000,
0.01570215782822780000,
0.01636209024170050000,
0.01703640867914670000,
0.01772593549749120000,
0.01842984628149530000,
0.01914814982886090000,
0.01988145867371130000,
0.02062890398582250000,
0.02139082000729940000,
0.02216739038276050000,
0.02295785988536490000,
0.02376282796984130000,
0.02458196525694420000,
0.02541469364786420000,
0.02626187536622570000,
0.02712268084763750000,
0.02799682962714990000,
0.02888510798855510000,
0.02978639071775980000,
0.03070093596244830000,
0.03162916227452560000,
0.03256977379897140000,
0.03352369359773240000,
0.03449052511261670000,
0.03546911580521190000,
0.03646104507071570000,
0.03746495715217230000,
0.03848011195762450000,
0.03950832508559470000,
0.04054766702274910000,
0.04159795210878940000,
0.04266055145239200000,
0.04373356025569980000,
0.04481724720337910000,
0.04591238822445110000,
0.04701715710708090000,
0.04813250843154150000,
0.04925814271146640000,
0.05039265858007830000,
0.05153752482854440000,
0.05269156739509920000,
0.05385379919716940000,
0.05502595683723440000,
0.05620614122190290000,
0.05739398298723640000,
0.05859090885782990000,
0.05979465555379200000,
0.06100598707332320000,
0.06222511219420360000,
0.06344991143268420000,
0.06468228446985730000,
0.06592093898213830000,
0.06716441904453470000,
0.06841507725737630000,
0.06967064658108250000,
0.07093019857889710000,
0.07219638021047560000,
0.07346605586640520000,
0.07473933326911790000,
0.07601805678222020000,
0.07729898025161000000,
0.07858342340749450000,
0.07987167069612440000,
0.08116090078220680000,
0.08245370230963710000,
0.08374853233628920000,
0.08504344881019380000,
0.08634142072196670000,
0.08763957244144220000,
0.08893700041385680000,
0.09023685004036040000,
0.09153527811269800000,
0.09283278177541530000,
0.09413091180643540000,
0.09542620818753040000,
0.09672082954190810000,
0.09801404545765440000,
0.09930328431921710000,
0.10059157052638000000,
0.10187658168111500000,
0.10315655725376600000,
0.10443527476815300000,
0.10570877895478100000,
0.10697680123321300000,
0.10824231779180900000,
0.10950091056749200000,
0.11075409200483100000,
0.11200281292229600000,
0.11324284196518500000,
0.11447808712318300000,
0.11570657544191100000,
0.11692461623093900000,
0.11813878095125400000,
0.11934331051821600000,
0.12053696200281300000,
0.12172570102238600000,
0.12290335646063700000,
0.12406998492892500000,
0.12522964615223800000,
0.12637654790548700000,
0.12751307577463900000,
0.12864011085088000000,
0.12975305040669500000,
0.13085604818609200000,
0.13194691509425400000,
0.13302270465267300000,
0.13408843651196900000,
0.13513975997071000000,
0.13617558872800900000,
0.13720059721283400000,
0.13820889671670500000,
0.13920234230328300000,
0.14018287461043100000,
0.14114423570393700000,
0.14209285939062300000,
0.14302560136039400000,
0.14393798990224600000,
0.14483806879875900000,
0.14571970165688200000,
0.14658040843863400000,
0.14742848726031000000,
0.14825575800259600000,
0.14906270304816900000,
0.14985504381827100000,
0.15062480219307300000,
0.15137551549362700000,
0.15210871261498000000,
0.15281779268967200000,
0.15350958341042700000,
0.15418088702498600000,
0.15482699640484600000,
0.15545703250327400000,
0.15606363865923100000,
0.15664504081355000000,
0.15721004523467400000,
0.15774934830898800000,
0.15826454512340800000,
0.15876099495204600000,
0.15922999143088000000,
0.15967773492773200000,
0.16010302387958200000,
0.16049974452252100000,
0.16087703098829600000,
0.16122910466903800000,
0.16155263747867100000,
0.16185688648348100000,
0.16213306350816600000,
0.16238183443297400000,
0.16261015446733100000,
0.16280802935673700000,
0.16298149717491000000,
0.16313204545778500000,
0.16324997228540000000,
0.16334766159726500000,
0.16341826386576000000,
0.16345523901514100000,
0.16347489021299700000,
0.16346377913774400000,
0.16342078534713500000,
0.16335886959488100000,
0.16326484351267100000,
0.16314172595479900000,
0.16299710338096300000,
0.16281922157465300000,
0.16261570843433000000,
0.16238747265927200000,
0.16212492457012300000,
0.16184067692018900000,
0.16152772759778200000,
0.16118120144939900000,
0.16081469500432500000,
0.16041640434090400000,
0.15998697029581400000,
0.15953704782015900000,
0.15905280538370400000,
0.15854269323225300000,
0.15800821923021700000,
0.15743839228163200000,
0.15684782204442700000,
0.15622932221778700000,
0.15557575727030900000,
0.15490445630997200000,
0.15420218144224000000,
0.15346691577922100000,
0.15271492962307900000,
0.15192903558745400000,
0.15111557990326000000,
0.15028251709980500000,
0.14941386485163800000,
0.14852460925819700000,
0.14761160443719400000,
0.14666116335331900000,
0.14569910259820400000,
0.14470726981107300000,
0.14368035913821700000,
0.14264406219434200000,
0.14157545454446200000,
0.14047733608022000000,
0.13936714425590800000,
0.13822380413504700000,
0.13705840763341800000,
0.13587623458445700000,
0.13466092958707500000,
0.13343119953680000000,
0.13217965948100800000,
0.13089661719319800000,
0.12960483359615800000,
0.12828707732803000000,
0.12694160421199600000,
0.12558969065898000000,
0.12420764356880800000,
0.12280862034488600000,
0.12139814966951200000,
0.11995611577392800000,
0.11850757578363400000,
0.11704215306888800000,
0.11554615596165600000,
0.11405112895786600000,
0.11253438647038200000,
0.11099127353207700000,
0.10945283168483300000,
0.10788885826293400000,
0.10630807956110200000,
0.10472837715736300000,
0.10312137422358100000,
0.10151052160666900000,
0.09989424771020730000,
0.09824833777688550000,
0.09661458275912110000,
0.09496663642024660000,
0.09329185070729760000,
0.09163767945950980000,
0.08996278939922200000,
0.08827192233316990000,
0.08659751127821650000,
0.08490202988147430000,
0.08320383621472270000,
0.08151486006120320000,
0.07980455651875050000,
0.07810715716205860000,
0.07640976198761880000,
0.07469340571073860000,
0.07300149862128520000,
0.07130239721175660000,
0.06959185543396300000,
0.06790996183488820000,
0.06621192027475380000,
0.06452233453408260000,
0.06285326182741120000,
0.06116377253543050000,
0.05950495935741000000,
0.05785520384198320000,
0.05618606089122140000,
0.05456267502158260000,
0.05293786465021880000,
0.05130179520674800000,
0.04971862264348100000,
0.04812379050466150000,
0.04653726777200880000,
0.04499569337889260000,
0.04343555507245920000,
0.04191230040370010000,
0.04041836903878840000,
0.03889715236115790000,
0.03745030471223890000,
0.03600872586292960000,
0.03454546437677380000,
0.03317248383283680000,
0.03178943862554570000,
0.03040846598878690000,
0.02910508822880960000,
0.02778491707733190000,
0.02649883458465700000,
0.02526903649007080000,
0.02401768345615770000,
0.02283766725310090000,
0.02168534391019460000,
0.02051482156561540000,
0.01944174787688840000,
0.01837275828262420000,
0.01730197315749950000,
0.01633393868407420000,
0.01534212125401040000,
0.01440090727908550000,
0.01352921602269550000,
0.01262361421948320000,
0.01181377104415380000,
0.01103839236351520000,
0.01023102348233120000,
0.00954566075804465000,
0.00886267945350381000,
0.00816929656619075000,
0.00759530044170959000,
0.00699772617182397000,
0.00643346227431199000,
0.00595158982967199000,
0.00543207835377655000,
0.00499616244929557000,
0.00459593988046467000,
0.00414942388751984000,
0.00382738027715532000,
0.00349801004821391000,
0.00314882950523732000,
0.00289473616639720000,
0.00262665752329835000,
0.00236589251516932000,
0.00216498915715933000,
0.00194779888001381000,
0.00176031892265548000,
0.00160161233568381000,
0.00142951496926075000,
0.00129711740126127000,
0.00117269670990650000,
0.00104000936321184000,
0.00094496786274736500,
0.00084747210882412900,
0.00075079490533337100,
0.00068062655441105000,
0.00060387363940618700,
0.00053835464649878100,
0.00048377497691693200,
0.00042478738009255000,
0.00037951307089179300,
0.00033828855355236800,
0.00029473137942714400,
0.00026323924900655500,
0.00023195052792712900,
0.00020093973956157600,
0.00017911863637681000,
0.00015470798791679500,
0.00013458192403953900,
0.00011890156259376100,
0.00010067762918919400,
0.00008847873303480690,
0.00007702878835965260,
0.00006434897490027950,
0.00005621762269842750,
0.00004836193143987330,
0.00004011353621237130,
0.00003470290343903300,
0.00002914883689138590,
0.00002402119970108930,
0.00002035413155423240,
0.00001658261242165370,
0.00001371579300650630,
0.00001147463960810510,
0.00000879243810874815,
0.00000732433394172841,
0.00000593452662368684,
0.00000483443481047313,
0.00000385720124558409,
0.00000287056272337702,
0.00000231651980810026,
0.00000173874539481987,
0.00000113272110508464,
0.00000117000390591950
}   ;
   
     for(G4int i = 0; i< 401 ; i++) {
      norma += ElectronEnergySp[i];
     }   
        
   for (G4int i = 0; i < 401; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;     
}

   for(G4int i = 0; i < 401; i++) {
     sum += ElectronEnergySp[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < 401; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
G4cout<<"ELECTRON ENERGY = "<<deltaX*y/x + EnergyBin[i]<<G4endl;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
   
   return 0;

}
//-------------------------------------------------------------------------
//    Energy Shooter for positrons from N13 beta+ decay
//-------------------------------------------------------------------------

G4double BxGeneratorNeutrinoC13::ShootEnergyPositronN13() {

   G4double sum = 0.;
   G4double norma = 0.;
   G4double Probability [401];
   G4double EnergyBin [401];
   G4double E0 = 1.201; //end Total Energy of the spectre
   G4double scale = (E0)/400;
   G4double xE;
   G4double deltaX, x, y;   

//Calculated by A. Ianni   
   G4double PositronEnergySp [401] =    
 {   
0.00279407988858406000,
0.03890662907939410000,
0.07840052335633310000,
0.11479749176822200000,
0.14822312794918400000,
0.17919735135988300000,
0.20816356450024200000,
0.23546265781212400000,
0.26135503028016400000,
0.28604222344721700000,
0.30968290257557100000,
0.33240408861661400000,
0.35430904174865500000,
0.37548286554727000000,
0.39599655836372300000,
0.41590999365387500000,
0.43527414946213100000,
0.45413280250271900000,
0.47252383407524700000,
0.49048025011051100000,
0.50803098758289400000,
0.52520155909630500000,
0.54201457334563100000,
0.55849015926768600000,
0.57464631466312600000,
0.59049919500076700000,
0.60606335441351600000,
0.62135194815985100000,
0.63637690378022200000,
0.65114906663498600000,
0.66567832433326100000,
0.67997371365700900000,
0.69404351288195500000,
0.70789532184733600000,
0.72153613169306500000,
0.73497238583900200000,
0.74821003350610000000,
0.76125457685797100000,
0.77411111266257900000,
0.78678436922792600000,
0.79927873924676000000,
0.81159830908690700000,
0.82374688498330000000,
0.83572801652005600000,
0.84754501773525500000,
0.85920098613389200000,
0.87069881985503800000,
0.88204123320580000000,
0.89323077074666300000,
0.90426982008859000000,
0.91516062354170200000,
0.92590528873804700000,
0.93650579833582300000,
0.94696401889931800000,
0.95728170903774700000,
0.96746052687674500000,
0.97750203692710500000,
0.98740771640899100000,
0.99717896108265300000,
1.00681709063147000000,
1.01632335363830000000,
1.02569893219149000000,
1.03494494615347000000,
1.04406245712109000000,
1.05305247210444000000,
1.06191594694775000000,
1.07065378951384000000,
1.07926686265177000000,
1.08775598696487000000,
1.09612194339551000000,
1.10436547564066000000,
1.11248729241180000000,
1.12048806955058000000,
1.12836845201173000000,
1.13612905572265000000,
1.14377046932916000000,
1.15129325583535000000,
1.15869795414543000000,
1.16598508051442000000,
1.17315512991405000000,
1.18020857731981000000,
1.18714587892456000000,
1.19396747328363000000,
1.20067378239602000000,
1.20726521272592000000,
1.21374215616856000000,
1.22010499096372000000,
1.22635408256065000000,
1.23248978443711000000,
1.23851243887569000000,
1.24442237769975000000,
1.25021992297189000000,
1.25590538765671000000,
1.26147907625055000000,
1.26694128537957000000,
1.27229230436878000000,
1.27753241578289000000,
1.28266189594132000000,
1.28768101540836000000,
1.29259003946023000000,
1.29738922853000000000,
1.30207883863209000000,
1.30665912176670000000,
1.31113032630633000000,
1.31549269736419000000,
1.31974647714648000000,
1.32389190528874000000,
1.32792921917738000000,
1.33185865425720000000,
1.33568044432545000000,
1.33939482181319000000,
1.34300201805485000000,
1.34650226354599000000,
1.34989578819058000000,
1.35318282153774000000,
1.35636359300875000000,
1.35943833211474000000,
1.36240726866547000000,
1.36527063296980000000,
1.36802865602784000000,
1.37068156971568000000,
1.37322960696281000000,
1.37567300192247000000,
1.37801199013532000000,
1.38024680868696000000,
1.38237769635930000000,
1.38440489377610000000,
1.38632864354310000000,
1.38814919038277000000,
1.38986678126422000000,
1.39148166552807000000,
1.39299409500678000000,
1.39440432414078000000,
1.39571261009006000000,
1.39691921284194000000,
1.39802439531491000000,
1.39902842345885000000,
1.39993156635152000000,
1.40073409629188000000,
1.40143628888996000000,
1.40203842315372000000,
1.40254078157296000000,
1.40294365020008000000,
1.40324731872851000000,
1.40345208056807000000,
1.40355823291814000000,
1.40356607683802000000,
1.40347591731538000000,
1.40328806333206000000,
1.40300282792787000000,
1.40262052826241000000,
1.40214148567460000000,
1.40156602574061000000,
1.40089447832969000000,
1.40012717765839000000,
1.39926446234277000000,
1.39830667544946000000,
1.39725416454452000000,
1.39610728174132000000,
1.39486638374650000000,
1.39353183190484000000,
1.39210399224241000000,
1.39058323550889000000,
1.38896993721804000000,
1.38726447768737000000,
1.38546724207639000000,
1.38357862042383000000,
1.38159900768371000000,
1.37952880376034000000,
1.37736841354225000000,
1.37511824693529000000,
1.37277871889446000000,
1.37035024945514000000,
1.36783326376325000000,
1.36522819210454000000,
1.36253546993306000000,
1.35975553789892000000,
1.35688884187503000000,
1.35393583298337000000,
1.35089696762033000000,
1.34777270748141000000,
1.34456351958515000000,
1.34126987629659000000,
1.33789225534994000000,
1.33443113987053000000,
1.33088701839640000000,
1.32726038489921000000,
1.32355173880447000000,
1.31976158501138000000,
1.31589043391204000000,
1.31193880141027000000,
1.30790720893972000000,
1.30379618348167000000,
1.29960625758242000000,
1.29533796936997000000,
1.29099186257055000000,
1.28656848652445000000,
1.28206839620172000000,
1.27749215221719000000,
1.27284032084542000000,
1.26811347403483000000,
1.26331218942200000000,
1.25843705034515000000,
1.25348864585748000000,
1.24846757074030000000,
1.24337442551553000000,
1.23820981645811000000,
1.23297435560808000000,
1.22766866078225000000,
1.22229335558574000000,
1.21684906942305000000,
1.21133643750907000000,
1.20575610087962000000,
1.20010870640191000000,
1.19439490678463000000,
1.18861536058776000000,
1.18277073223238000000,
1.17686169200991000000,
1.17088891609147000000,
1.16485308653667000000,
1.15875489130254000000,
1.15259502425201000000,
1.14637418516223000000,
1.14009307973284000000,
1.13375241959379000000,
1.12735292231323000000,
1.12089531140513000000,
1.11438031633666000000,
1.10780867253547000000,
1.10118112139671000000,
1.09449841029010000000,
1.08776129256660000000,
1.08097052756507000000,
1.07412688061876000000,
1.06723112306160000000,
1.06028403223435000000,
1.05328639149079000000,
1.04623899020345000000,
1.03914262376956000000,
1.03199809361653000000,
1.02480620720768000000,
1.01756777804745000000,
1.01028362568683000000,
1.00295457572847000000,
0.99558145983178000000,
0.98816511571779800000,
0.98070638717416600000,
0.97320612405975900000,
0.96566518230937100000,
0.95808442393825100000,
0.95046471704655600000,
0.94280693582366400000,
0.93511196055245600000,
0.92738067761350300000,
0.91961397948910300000,
0.91181276476727500000,
0.90397793814571500000,
0.89611041043557500000,
0.88821109856523000000,
0.88028092558393900000,
0.87232082066545900000,
0.86433171911154800000,
0.85631456235541100000,
0.84827029796510100000,
0.84019987964675600000,
0.83210426724795600000,
0.82398442676079700000,
0.81584133032502400000,
0.80767595623110100000,
0.79948928892318500000,
0.79128231900208800000,
0.78305604322800500000,
0.77481146452352900000,
0.76654959197624800000,
0.75827144084149300000,
0.74997803254500200000,
0.74167039468546000000,
0.73334956103713100000,
0.72501657155223900000,
0.71667247236349900000,
0.70831831578645000000,
0.69995516032184400000,
0.69158407065791900000,
0.68320611767269400000,
0.67482237843612000000,
0.66643393621234000000,
0.65804188046170400000,
0.64964730684294800000,
0.64125131721523800000,
0.63285501964010600000,
0.62445952838348800000,
0.61606596391761800000,
0.60767545292297700000,
0.59928912829007700000,
0.59090812912134300000,
0.58253360073286600000,
0.57416669465619500000,
0.56580856864003500000,
0.55746038665194800000,
0.54912331887999800000,
0.54079854173440200000,
0.53248723784911900000,
0.52419059608340700000,
0.51590981152335900000,
0.50764608548346500000,
0.49940062550800700000,
0.49117464537260200000,
0.48296936508557600000,
0.47478601088937600000,
0.46662581526196300000,
0.45849001691815100000,
0.45037986081093000000,
0.44229659813278700000,
0.43424148631697500000,
0.42621578903875100000,
0.41822077621664000000,
0.41025772401363700000,
0.40232791483837400000,
0.39443263734632400000,
0.38657318644093000000,
0.37875086327471800000,
0.37096697525041300000,
0.36322283602205500000,
0.35551976549601100000,
0.34785908983207200000,
0.34024214144446500000,
0.33267025900284500000,
0.32514478743332300000,
0.31766707791941200000,
0.31023848790301500000,
0.30286038108533700000,
0.29553412742783100000,
0.28826110315310500000,
0.28104269074581000000,
0.27388027895351200000,
0.26677526278756200000,
0.25972904352396000000,
0.25274302870414500000,
0.24581863213586200000,
0.23895727389393600000,
0.23216038032106300000,
0.22542938402859800000,
0.21876572389731600000,
0.21217084507815700000,
0.20564619899297100000,
0.19919324333522400000,
0.19281344207073800000,
0.18650826543836900000,
0.18027918995070200000,
0.17412769839473500000,
0.16805527983252800000,
0.16206342960187000000,
0.15615364931691200000,
0.15032744686881000000,
0.14458633642632900000,
0.13893183843647300000,
0.13336547962507000000,
0.12788879299736200000,
0.12250331783859700000,
0.11721059971459100000,
0.11201219047228700000,
0.10690964824030800000,
0.10190453742950100000,
0.09699842873346640000,
0.09219289912908720000,
0.08748953187703870000,
0.08288991652229660000,
0.07839564889463670000,
0.07400833110912200000,
0.06972957156658910000,
0.06556098495411820000,
0.06150419224549590000,
0.05756082070168320000,
0.05373250387125370000,
0.05002088159084810000,
0.04642759998560080000,
0.04295431146957420000,
0.03960267474617740000,
0.03637435480858110000,
0.03327102294012690000,
0.03029435671472620000,
0.02744603999725470000,
0.02472776294394230000,
0.02214122200275200000,
0.01968811991375700000,
0.01737016570950900000,
0.01518907471540110000,
0.01314656855002540000,
0.01124437512552380000,
0.00948422864793394000,
0.00786786961752842000,
0.00639704482914938000,
0.00507350737253756000,
0.00389901663265584000,
0.00287533829000654000,
0.00200424432094549000,
0.00128751299798913000,
0.00072692889011809000,
0.00032428286307435500,
0.00008137207965495030,
0.00000000000000000000
}
;


   
     for(G4int i = 0; i< 401 ; i++) {
      norma += PositronEnergySp[i];
     }      

   for (G4int i = 0; i < 401; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;
}

   for(G4int i = 0; i < 401; i++) {
     sum += PositronEnergySp[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < 401; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
                     
G4cout<<"POSITRON ENERGY = "<<deltaX*y/x + EnergyBin[i]<<G4endl;
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }
   
   return 0;
}


G4double BxGeneratorNeutrinoC13::ShootEnergyElectronFlat() {

   G4double sum = 0.;
   G4double norma = 0.;
   G4double Probability [401];
   G4double EnergyBin [401];
   G4double E0 = 15.0; //end Total Energy of the spectre
   G4double scale = E0/400;
   G4double xE;

  G4double deltaX, x, y;   

  G4double ElectronEnergySp [401] ;
    

     for(G4int i = 0; i< 401 ; i++) {
       ElectronEnergySp[i] = 1 ;
     }      

      ElectronEnergySp[0] = 0 ;
 
     for(G4int i = 0; i< 401 ; i++) {
      norma += ElectronEnergySp[i];
     }      

   for (G4int i = 0; i < 401; i++) {
     xE = scale*G4float(i);
     EnergyBin[i] = xE;
}

   for(G4int i = 0; i < 401; i++) {
     sum += ElectronEnergySp[i]/norma;
     Probability[i] = sum;
   }

   G4double val = G4UniformRand();

   for(G4int i = 0; i < 401; i++) {
     if(Probability[i] >= val) {
       if(i == 0) return EnergyBin[0];
       deltaX = val - Probability[i] ;
       y = EnergyBin[i] - EnergyBin[i-1] ;
       x = Probability[i] - Probability[i-1] ;
                     
       return (deltaX*y/x + EnergyBin[i]) ;
     }
   }   
   return 0;
}
