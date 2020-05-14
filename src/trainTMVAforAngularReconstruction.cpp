/* \file trainTMVAforAngularReconstruction
   \brief use TMVA methods for modified disp method

   this code is used for training for angular and energy reconstruction


*/

#include "TChain.h"
#include "TCut.h"
#include "TFile.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Config.h"
#ifdef ROOT6
#include "TMVA/DataLoader.h"
#endif
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Ctelconfig.h"
#include "Cshowerpars.h"
#include "Ctpars.h"
#include "VCameraRead.h"
#include "VEmissionHeightCalculator.h"
#include "VDetectorTree.h"
#include "VGlobalRunParameter.h"
#include "VUtilities.h"

using namespace std;

/////////////////////////////////////////////////////
// one tree per telescope type
map< ULong64_t, TTree* > fMapOfTrainingTree;
/////////////////////////////////////////////////////

/*

    train TVMA method and write results into the corresponding directory

    one MVA per telescope type

    Allowed target BDTs:

    BDTDisp
    BDTDispError
    BDTDispEnergy
    BDTDispCore


*/
bool trainTMVA( string iOutputDir, float iTrainTest,
                ULong64_t iTelType, TTree* iDataTree,
                string iTargetBDT, string iTMVAOptions,
                string iQualityCut )
{
    cout << endl;
    cout << "Starting " << iTargetBDT;
    cout << " training for telescope type " << iTelType << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << endl;
    
    if( !iDataTree )
    {
        cout << "Error: data tree for telescope type " << iTelType << " does not exist" << endl;
        return false;
    }
    char htitle[6000];
    
    // Determine the number of training and test events
    unsigned int ntrain   = 0 ;
    unsigned int ntest    = 0 ;
    unsigned int nentries = iDataTree->GetEntries() ;
    cout << endl;
    ntrain = floor( nentries * iTrainTest ) ;
    ntest  = nentries - ntrain ;
    if( ntrain <= 100 )
    {
        cout << endl;
        cout << "Error, 4th argument train/test fraction is so small that only " << ntrain;
        cout << " events were selected for training, while TMVA usually needs thousands of training events to work properly.";
        cout << "Try increasing the 4th argument... (you only have " << nentries;
        cout << " total events to designate for either training or testing...)" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    if( ntest <= 100 )
    {
        cout << endl;
        cout << "Error, 4th argument train/test fraction is so large that only " << ntest;
        cout << " events were selected for testing, while TMVA usually needs thousands of testing events to work properly.";
        cout << " Try decreasing the 4th argument... (you only have " << nentries;
        cout << " total events to designate for either training or testing...)" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    ntrain *= 0.8;
    ntest *= 0.8;
    cout << "\tnumber of training events: " << ntrain << endl;
    cout << "\tnumber of test events    : " << ntest  << endl;
    cout << "\tfraction of test events  : " << iTrainTest << endl << endl;
    std::ostringstream train_and_test_conditions ;
    train_and_test_conditions << "nTrain_Regression=" << ntrain << ":" ;
    train_and_test_conditions << "nTest_Regression="  << ntest  << ":" ;
    train_and_test_conditions << "SplitMode=Random"             << ":" ;
    train_and_test_conditions << "NormMode=NumEvents"           << ":" ;
#ifdef ROOT6
    train_and_test_conditions << "V=True" << ":" ;
    train_and_test_conditions << "VerboseLevel=Info" << ":" ;
    train_and_test_conditions << "ScaleWithPreselEff=True";
#else
    train_and_test_conditions << "!V";
#endif
    cout << "Train and test condition: " << train_and_test_conditions.str() << endl;
    cout << endl;
    
    
    // output file name
    ostringstream iFileName;
    iFileName << iOutputDir << "/" << iTargetBDT << "_" << iTelType << ".tmva.root";
    TFile* i_tmva = new TFile( iFileName.str().c_str(), "RECREATE" );
    if( i_tmva->IsZombie() )
    {
        cout << "error while creating tmva root file: " << iFileName.str() << endl;
        return false;
    }
    
    // set output directory
    gSystem->mkdir( iOutputDir.c_str() );
    TString iOutputDirectory( iOutputDir.c_str() );
    gSystem->ExpandPathName( iOutputDirectory );
    ( TMVA::gConfig().GetIONames() ).fWeightFileDir = iOutputDirectory;
    
    
    // tmva regression
#ifdef ROOT6
    TMVA::Factory* factory = new TMVA::Factory( iTargetBDT.c_str(), i_tmva, 
                            "V:!DrawProgressBar:!Color:!Silent:AnalysisType=Regression:VerboseLevel=Debug:Correlations=True" );
    factory->SetVerbose( true );
#else
    TMVA::Factory* factory = new TMVA::Factory( iTargetBDT.c_str(), i_tmva, 
                            "V:!DrawProgressBar:!Color:!Silent" );
#endif
    
#ifdef ROOT6
    TMVA::DataLoader* dataloader = new TMVA::DataLoader( "" );
#else
    TMVA::Factory* dataloader = factory;
#endif
    
    // list of variables used by MVA method
    dataloader->AddVariable( "width" , 'F' );
    dataloader->AddVariable( "length", 'F' );
    dataloader->AddVariable( "wol",    'F' );
    dataloader->AddVariable( "size"  , 'F' );
    // hard coded ASTRI telescope type
    // (no time gradient is available)
    if( iTelType != 201511619 )
    {
        dataloader->AddVariable( "tgrad_x*tgrad_x", 'F' );
    }
    //No stereo variables used
    //dataloader->AddVariable( "cross" , 'F' );
    dataloader->AddVariable( "asym"  , 'F' );
    dataloader->AddVariable( "loss"  , 'F' );
    dataloader->AddVariable( "dist"  , 'F' );
    dataloader->AddVariable( "fui"  , 'F' );
    /*if( iTargetBDT == "BDTDispEnergy" )
    {
        dataloader->AddVariable( "EHeight", 'F' );
        dataloader->AddVariable( "Rcore", 'F' );
    }*/
    // spectators
    dataloader->AddSpectator( "MCe0", 'F' );
    dataloader->AddSpectator( "MCxoff", 'F' );
    dataloader->AddSpectator( "MCyoff", 'F' );
    dataloader->AddSpectator( "MCxcore", 'F' );
    dataloader->AddSpectator( "MCycore", 'F' );
    dataloader->AddSpectator( "MCrcore", 'F' );
    dataloader->AddSpectator( "NImages", 'F' );
    // train for energy reconstruction
    if( iTargetBDT == "BDTDispEnergy" )
    {
        // dispEnergy is log10(E)/log10(size)
        dataloader->AddTarget( "dispEnergy", 'F' );
    }
    // train for direction reconstruction
    else if( iTargetBDT == "BDTDisp" )
    {
        dataloader->AddTarget( "disp"  , 'F' );
    }
    // train for error on disp reconstruction
    else if( iTargetBDT == "BDTDispError" )
    {
        dataloader->AddTarget( "dispError", 'F', "dispError", 0., 10. );
    }
    // train for core reconstruction
    else if( iTargetBDT == "BDTDispCore" )
    {
        dataloader->AddTarget( "dispCore"  , 'F', "m", 0., 1.e5 );
    }
    else
    {
        cout << "Error: unknown target BDT: " << iTargetBDT << endl;
        exit( EXIT_FAILURE );
    }
    // add weights (optional)
    //    dataloader->SetWeightExpression( "MCe0*MCe0", "Regression" );
    
    // regression tree
    dataloader->AddRegressionTree( iDataTree, 1. );
    
    // quality cuts
    // (determined by plotting all variables with
    //  macro plot_dispBDT_inputVariables.C)
    //  **IMPORTANT**
    //  loss cut here must correspond to loss cut later in the analysis
    //  (otherwise large bias in the energy reconstruction)
    TCut fQualityCut = iQualityCut.c_str();
    cout << "Quality cuts applied: " << iQualityCut << endl;
    
    dataloader->PrepareTrainingAndTestTree( fQualityCut, train_and_test_conditions.str().c_str() ) ;
    
    ostringstream iMVAName;
    iMVAName << "BDT_" << iTelType;
    sprintf( htitle, "%s", iMVAName.str().c_str() );
    
    TString methodstr( iTMVAOptions.c_str() );
    
    cout << "Built MethodStringStream: " << iTMVAOptions << endl;
    cout << endl;
    TString methodTitle( htitle );
#ifdef ROOT6
    factory->BookMethod( dataloader, TMVA::Types::kBDT, methodTitle, methodstr ) ;
#else
    factory->BookMethod( TMVA::Types::kBDT, methodTitle, methodstr ) ;
#endif
    
    factory->TrainAllMethods();
    
    factory->TestAllMethods();
    
    factory->EvaluateAllMethods();
    
    factory->Delete();
    
    return true;
}

/*
 * read a ascii file with a list of evndisplay files
 *
 * return a vector with all the files
 *
*/
vector< string > fillInputFile_fromList( string iList )
{
    vector< string > inputfile;
    
    
    ifstream is;
    is.open( iList.c_str(), ifstream::in );
    if( !is )
    {
        cout << "fillInputFile_fromList() error reading list of input files: " << endl;
        cout << iList << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    cout << "fillInputFile_fromList() reading input file list: ";
    cout << iList << endl;
    string iLine;
    while( getline( is, iLine ) )
    {
        if( iLine.size() > 0 )
        {
            inputfile.push_back( iLine );
        }
    }
    is.close();
    
    cout << "total number of input files " << inputfile.size() << endl;
    
    return inputfile;
}

/******************************************************************************************

    reading training data used for the TMVA training

    (a previous training session produced these files
*/
bool readTrainingFile( string iTargetBDT, ULong64_t iTelType, const string iDataDirectory )
{
    fMapOfTrainingTree.clear();
    
    ostringstream iFileName;
    iFileName << iDataDirectory << "/" << iTargetBDT;
    if( iTelType != 0 )
    {
        iFileName << "_" << iTelType;
    }
    iFileName << ".root";
    TFile* iO = new TFile( iFileName.str().c_str() );
    if( iO->IsZombie() )
    {
        cout << "error reading training trees from file ";
        cout << iFileName.str() << endl;
        return false;
    }
    ostringstream iTreeName;
    iTreeName << "dispTree_" << iTelType;
    if( iO->Get( iTreeName.str().c_str() ) )
    {
        fMapOfTrainingTree[iTelType] = ( TTree* )iO->Get( iTreeName.str().c_str() );
    }
    
    return true;
}

/*
 * read list of valid telescopes from a typical telescope (array) list
 *
 */
vector< bool > readArrayList( unsigned int i_ntel, string iArrayList, vector< unsigned int > iHyperArrayID )
{
    vector< bool > iTelList( i_ntel, true );
    
    if( iArrayList.size() > 0 )
    {
        // switch all telescopes off
        iTelList.assign( iTelList.size(), false );
        
        // open telescope list file
        ifstream is;
        is.open( iArrayList.c_str(), ifstream::in );
        if( !is )
        {
            cout << "readArrayList() error reading list of arrays from " << endl;
            cout << iArrayList << endl;
            cout << "exiting..." << endl;
            exit( EXIT_FAILURE );
        }
        cout << "reading list of telescope from " << iArrayList << endl;
        string iLine;
        while( getline( is, iLine ) )
        {
            if( iLine.size() > 0 )
            {
                istringstream is_stream( iLine );
                unsigned int T = ( unsigned int )atoi( is_stream.str().c_str() );
                // find position of telescope in hyper array list
                for( unsigned int i = 0; i < iHyperArrayID.size(); i++ )
                {
                    if( T == iHyperArrayID[i] )
                    {
                        if( i < iTelList.size() )
                        {
                            iTelList[i] = true;
                        }
                    }
                }
            }
        }
    }
    return iTelList;
}



/******************************************************************************************

   write the training file used for the TMVA training
   (simply a tree with all the necessary variables; one tree per telescope type)

   output as a root file (might be temporary, steer with scripts)

*/
bool writeTrainingFile( const string iInputFile, ULong64_t iTelType,
                        unsigned int iRecID, string iArrayList )
{
    ////////////////////////////
    // read list of input files
    vector< string > iInputFileList = fillInputFile_fromList( iInputFile );
    if( iInputFileList.size() == 0 )
    {
        cout << "writeTrainingFile error: input file list is empty" << endl;
        return false;
    }
    // (assume in the following that iInputFileList has a reasonable size)
    
    /////////////////////////////
    // telescope configuration
    TChain i_telChain( "telconfig" );
    i_telChain.Add( iInputFileList[0].c_str(), 0 );
    
    Ctelconfig i_tel( &i_telChain );
    i_tel.GetEntry( 0 );
    unsigned int i_ntel = i_tel.NTel;
    
    vector< unsigned int > iHyperArrayID;
    vector< float > iFOV_tel;
    
    // get list of telescopes - hyperarray values
    for( int t = 0; t < i_tel.fChain->GetEntries(); t++ )
    {
        i_tel.GetEntry( t );
        
        iHyperArrayID.push_back( i_tel.TelID_hyperArray );
        iFOV_tel.push_back( i_tel.FOV );
        cout << "\t FOV for telescope " << iHyperArrayID.back() << ": " << iFOV_tel.back() << endl;
    }
    
    // read list of telescope from usual array lists
    vector< bool > fUseTelescope = readArrayList( i_ntel, iArrayList, iHyperArrayID );
    if( fUseTelescope.size() != i_ntel )
    {
        cout << "Error in telescope list size: " << fUseTelescope.size() << ", " << i_ntel << endl;
        cout << "Exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    
    //
    // vector with telescope position
    // (includes all telescopes, even those
    // of other types)
    vector< float > fTelX;
    vector< float > fTelY;
    vector< float > fTelZ;
    vector< ULong64_t > fTelType;
    unsigned int f_ntelType = 0;
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        i_tel.GetEntry( i );
        
        fTelX.push_back( i_tel.TelX );
        fTelY.push_back( i_tel.TelY );
        fTelZ.push_back( i_tel.TelZ );
        fTelType.push_back( i_tel.TelType );
        
        if( i < fUseTelescope.size() && !fUseTelescope[i] )
        {
            continue;
        }
        if( i_tel.TelType == iTelType 
          || iTelType == 0 )
        {
            f_ntelType++;
        }
    }
    
    ///////////////////////////////////////////////////
    // definition of training trees (one per telescope type)
    
    int runNumber = -1;
    int eventNumber = -1;
    unsigned int tel = 0;
    float cen_x = -1.;
    float cen_y = -1.;
    float sinphi = -1.;
    float cosphi = -1.;
    float size = -1.;    // actually log10(size)
    float ntubes = -1.;
    float loss = -1.;
    float asym = -1.;
    float width = -1.;
    float length = -1.;
    float wol = -1.;    // width over length
    float MCe0 = -1.;
    float MCxoff = -1.;
    float MCyoff = -1.;
    float MCxcore = -1.;
    float MCycore = -1.;
    float MCrcore = -1.;
    float Xcore = -1.;
    float Ycore = -1.;
    float Rcore = -1.;
    float Xoff = -1.;
    float Yoff = -1.;
    float LTrig = -1.;
    float MCaz = -1.;
    float MCze = -1.;
    float disp = -1.;
    float dispError = -1.;
    float NImages = -1.;
    float cross = -1.;
    float dispPhi = -1.;
    float dispEnergy = -1.;
    float dispCore = -1.;
    float dist = -1.;
    float fui = -1.;
    float tgrad_x = -1.;
    float meanPedvar_Image = -1.;
    float ze = -1.;
    float az = -1.;
    float EmissionHeight = -1.;
    
    fMapOfTrainingTree.clear();
    cout << "total number of telescopes: " << i_ntel;
    cout << " (selected " << f_ntelType << ")" << endl;
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        i_tel.GetEntry( i );
        
        // select telescope type
        if( iTelType != 0 && i_tel.TelType != iTelType )
        {
            continue;
        }
        // check if telescope is in telescope list
        if( i < fUseTelescope.size() && !fUseTelescope[i] )
        {
            continue;
        }
        
        if( fMapOfTrainingTree.find( i_tel.TelType ) == fMapOfTrainingTree.end() )
        {
            ostringstream iTreeName;
            iTreeName << "dispTree_" << i_tel.TelType;
            ostringstream iTreeTitle;
            iTreeTitle << "training tree for modified disp method (telescope type " << i_tel.TelType << ")";
            fMapOfTrainingTree[i_tel.TelType] = new TTree( iTreeName.str().c_str(), iTreeTitle.str().c_str() );
            
            fMapOfTrainingTree[i_tel.TelType]->Branch( "runNumber"  , &runNumber  , "runNumber/I" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "eventNumber", &eventNumber, "eventNumber/I" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "tel",         &tel        , "tel/i" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cen_x"      , &cen_x      , "cen_x/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cen_y"      , &cen_y      , "cen_y/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "sinphi"     , &sinphi     , "sinphi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cosphi"     , &cosphi     , "cosphi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "size"       , &size       , "size/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "ntubes"     , &ntubes     , "ntubes/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "loss"       , &loss       , "loss/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "asym"       , &asym       , "asym/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "width"      , &width      , "width/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "length"     , &length     , "length/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "wol"        , &wol        , "wol/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dist"       , &dist       , "dist/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "fui"        , &fui        , "fui/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "tgrad_x"    , &tgrad_x    , "tgrad_x/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "meanPedvar_Image"         , &meanPedvar_Image, "meanPedvar_Image/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCe0"       , &MCe0       , "MCe0/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCxoff"     , &MCxoff     , "MCxoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCyoff"     , &MCyoff     , "MCyoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCxcore"    , &MCxcore    , "MCxcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCycore"    , &MCycore    , "MCycore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCrcore"    , &MCrcore    , "MCrcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Xcore"      , &Xcore      , "Xcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Ycore"      , &Ycore      , "Ycore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Rcore"      , &Rcore      , "Rcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Xoff"       , &Xoff       , "Xoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Yoff"       , &Yoff       , "Yoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "LTrig"      , &LTrig      , "LTrig/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "NImages"    , &NImages    , "NImages/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "EHeight"    , &EmissionHeight, "EHeight/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCaz"       , &MCaz       , "MCaz/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCze"       , &MCze       , "MCze/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Ze"         , &ze         , "Ze/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Az"         , &az         , "Az/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "disp"       , &disp       , "disp/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispError"  , &dispError  , "dispError/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cross"      , &cross      , "cross/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispPhi"    , &dispPhi    , "dispPhi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispEnergy" , &dispEnergy , "dispEnergy/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispCore"   , &dispCore   , "dispCore/F" );
            
        }
    }
    
    ////////////////////////////////////////////
    // filling of training trees;
    cout << "filling training trees for " << fMapOfTrainingTree.size() << " telescope type(s)" << endl;
    cout << "\t found " << f_ntelType << " telescopes of telescope type " << iTelType << endl;
    
    // get showerpars tree
    TChain i_showerparsTree( "showerpars" );
    for( unsigned int f = 0; f < iInputFileList.size(); f++ )
    {
        i_showerparsTree.Add( iInputFileList[f].c_str(), 0 );
    }
    Cshowerpars i_showerpars( &i_showerparsTree, true, true );
    
    // get all tpars tree
    vector< TChain* > i_tparsTree;
    vector< Ctpars* > i_tpars;
    for( unsigned int i = 0; i < fTelType.size(); i++ )
    {
        if( iTelType == 0 || iTelType == fTelType[i] )
        {
            ostringstream iTreeName;
            iTreeName << "Tel_" << i + 1 << "/tpars";
            i_tparsTree.push_back( new TChain( iTreeName.str().c_str() ) );
            for( unsigned int f = 0; f < iInputFileList.size(); f++ )
            {
                i_tparsTree.back()->Add( iInputFileList[f].c_str(), 0 );
            }
            i_tpars.push_back( new Ctpars( i_tparsTree.back(), true, true ) );
            cout << "\t found tree " << iTreeName.str();
            cout << " (teltype " << fTelType[i] << ")";
            cout << ", entries: ";
            cout << i_tpars.back()->fChain->GetEntries() << endl;
        }
        else
        {
            i_tpars.push_back( 0 );
            cout << "\t ignore tree for telescope type " << fTelType[i] << endl;
        }
    }
    
    // temporary variables for emission height calculation
    VEmissionHeightCalculator* fEmissionHeightCalculator = new VEmissionHeightCalculator();
    fEmissionHeightCalculator->setTelescopePositions( fTelX, fTelY, fTelZ );
    double fEM_cen_x[fTelType.size()];
    double fEM_cen_y[fTelType.size()];
    double fEM_size[fTelType.size()];
    
    /////////////////////////////////////////////////
    // loop over all events in trees
    int nentries = i_showerpars.fChain->GetEntries();
    cout << "Loop over " << nentries << " entries in source files" << endl;
    
    for( int n = 0; n < nentries; n++ )
    {
    	// read events from event trees
        i_showerpars.GetEntry( n );
        
        // check recid
       /* if( iRecID >= i_showerpars.NMethods )
        {
            cout << "Error: invalid reconstruction ID.";
            cout << " Maximum allowed value is " << i_showerpars.NMethods << endl;
            return false;
        }
        
        // require:
        // - reconstructed event
        // - at least two telescopes
        if( i_showerpars.Chi2[iRecID] < 0.
                ||  i_showerpars.NImages[iRecID] < 2 )
        {
            continue;
        }*/
        
        // check if there are image of this particle teltype
        /*int i_nteltypecounter = 0;
        for( unsigned int i = 0; i < fTelType.size(); i++ )
        {
            if( (fTelType[i] == iTelType || iTelType == 0)
                    && ( int )i_showerpars.ImgSel_list[iRecID][i] > 0 )
            {
                i_nteltypecounter++;
            }
        }
        if( i_nteltypecounter == 0 )
        {
            continue;
        }*/
        
        /////////////////////////////////////////////////////////
        // calculate emission height - unused in this analysis
        /*for( unsigned int i = 0; i < i_tpars.size(); i++ )
        {
            fEM_size[i] = -1.;
            fEM_cen_x[i] = 0.;
            fEM_cen_y[i] = 0.;
            
            if( ( int )i_showerpars.ImgSel_list[iRecID][i] < 1
                    && i_showerpars.NImages[iRecID] > 1 )
            {
                continue;
            }
            // note: this is different to what happens in the analysis:
            // here the emissionheight is calculated only from the
            // telescopes of the given telescope type
            if( !i_tpars[i] )
            {
                continue;
            }
            
            i_tpars[i]->GetEntry( n );
            
            if( i_tpars[i]->size > 0. )
            {
                fEM_size[i] = i_tpars[i]->size;
                fEM_cen_x[i] = i_tpars[i]->cen_x;
                fEM_cen_y[i] = i_tpars[i]->cen_y;
            }
        }
        EmissionHeight = fEmissionHeightCalculator->getEmissionHeight( fEM_cen_x, fEM_cen_y, fEM_size,
                         i_showerpars.ArrayPointing_Azimuth,
                         i_showerpars.ArrayPointing_Elevation );
*/
        //////////////////////////////////////
        // loop over all telescopes
        for( unsigned int i = 0; i < i_tpars.size(); i++ )
        {
            // check if telescope was reconstructed
            /*if( ( int )i_showerpars.ImgSel_list[iRecID][i] < 1
                    || i_showerpars.NImages[iRecID] < 2 )
            {
                continue;
            }
            // check if telescope is of valid telescope type
            if( (fTelType[i] != iTelType && iTelType != 0 ) || !i_tpars[i] )
            {
                continue;
            }
            // check if event is not completely out of the FOV
            // (use 20% x size of the camera)
            if( i < iFOV_tel.size()
                    && sqrt( i_showerpars.MCxoff * i_showerpars.MCxoff + i_showerpars.MCyoff * i_showerpars.MCyoff ) > iFOV_tel[i] * 0.5 * 1.2 )
            {
                continue;
            }*/
            i_tpars[i]->GetEntry( n );
            
            /////////////////////////////////
            // quality cuts
            if( i_tpars[i]->size <= 0 )
            {
                continue;
            }
            
            runNumber   = i_showerpars.runNumber;
            eventNumber = i_showerpars.eventNumber;
            tel         = i + 1;
            cen_x       = i_tpars[i]->cen_x;
            cen_y       = i_tpars[i]->cen_y;
            sinphi      = i_tpars[i]->sinphi;
            cosphi      = i_tpars[i]->cosphi;
            size        = log10( i_tpars[i]->size );
            ntubes      = i_tpars[i]->ntubes;
            loss        = i_tpars[i]->loss;
            asym        = i_tpars[i]->asymmetry;
            width       = i_tpars[i]->width;
            length      = i_tpars[i]->length;
            if( length > 0. )
            {
                wol = width / length;
            }
            else
            {
                wol = 0.;
            }
            dist        = i_tpars[i]->dist;
            fui         = i_tpars[i]->fui;
            tgrad_x     = i_tpars[i]->tgrad_x;
            meanPedvar_Image = i_tpars[i]->meanPedvar_Image;
            ze          = 90. - i_showerpars.TelElevation[i];
            az          = i_showerpars.TelAzimuth[i];
            MCe0        = i_showerpars.MCe0;
            MCxoff      = i_showerpars.MCxoff;
            MCyoff      = i_showerpars.MCyoff;
            MCxcore     = i_showerpars.MCxcore;
            MCycore     = i_showerpars.MCycore;
            Xoff        = i_showerpars.Xoff[iRecID];
            Yoff        = i_showerpars.Yoff[iRecID];
            Xcore       = i_showerpars.Xcore[iRecID];
            Ycore       = i_showerpars.Ycore[iRecID];
            LTrig       = i_showerpars.LTrig;
            NImages     = i_showerpars.NImages[iRecID];
            MCze        = i_showerpars.MCze;
            MCaz        = i_showerpars.MCaz;
            
            Rcore       = VUtilities::line_point_distance( Ycore,   -1.*Xcore,   0., ze, az, fTelY[i], -1.*fTelX[i], fTelZ[i] );
            MCrcore     = VUtilities::line_point_distance( MCycore, -1.*MCxcore, 0., MCze, MCaz, fTelY[i], -1.*fTelX[i], fTelZ[i] );
            
            /*if( Rcore < 0. )
            {
                continue;
            }*/
            
            //////////////////////////////////////////////////////////////////////////////////////////////////
            // calculate disp (observe sign convention for MC in y direction for MCyoff and Yoff)
            disp  = sqrt( ( cen_y + MCyoff ) * ( cen_y + MCyoff ) + ( cen_x - MCxoff ) * ( cen_x - MCxoff ) );
            cross = sqrt( ( cen_y + Yoff ) * ( cen_y + Yoff ) + ( cen_x - Xoff ) * ( cen_x - Xoff ) );
            dispPhi = TMath::ATan2( sinphi, cosphi ) - TMath::ATan2( cen_y + MCyoff, cen_x - MCxoff );
            
            // disp error: the expected difference between true and
            //             reconstructed direction
            dispError = 0;
            float x1 = cen_x - disp * cosphi;
            float x2 = cen_x + disp * cosphi;
            float y1 = cen_y - disp * sinphi;
            float y2 = cen_y + disp * sinphi;
            if( sqrt( ( x1 - MCxoff ) * ( x1 - MCxoff ) + ( y1 + MCyoff ) * ( y1 + MCyoff ) )
                    < sqrt( ( x2 - MCxoff ) * ( x2 - MCxoff ) + ( y2 + MCyoff ) * ( y2 + MCyoff ) ) )
            {
                dispError = sqrt( ( x1 - MCxoff ) * ( x1 - MCxoff ) + ( y1 + MCyoff ) * ( y1 + MCyoff ) );
            }
            else
            {
                dispError = sqrt( ( x2 - MCxoff ) * ( x2 - MCxoff ) + ( y2 + MCyoff ) * ( y2 + MCyoff ) );
            }
            
            // training target in ratio to size
            dispEnergy = log10( i_showerpars.MCe0 ) / log10( i_tpars[i]->size );
            dispCore   = Rcore;
            
            if( fMapOfTrainingTree.find( fTelType[i] ) != fMapOfTrainingTree.end() )
            {
                fMapOfTrainingTree[fTelType[i]]->Fill();
            }
        }
    }
    // cleanup
    for( unsigned int i = 0; i < i_tpars.size(); i++ )
    {
        if( i_tpars[i] )
        {
            delete i_tpars[i];
        }
    }
    
    return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

int main( int argc, char* argv[] )
{
    // print version only
    if( argc == 2 )
    {
        string fCommandLine = argv[1];
        if( fCommandLine == "-v" || fCommandLine == "--version" )
        {
            VGlobalRunParameter fRunPara;
            cout << fRunPara.getEVNDISP_VERSION() << endl;
            exit( EXIT_SUCCESS );
        }
    }
    cout << endl;
    
    // print help text
    if( argc < 6 )
    {
        cout << "./trainTMVAforAngularReconstruction <list of input eventdisplay files (MC)> <output directory>";
        cout << " <train vs test fraction> <RecID> <telescope type> [train for angular / energy / core reconstruction]";
        cout << " [MVA options] [array layout file] [directory with training trees] [quality cut]" << endl;
        cout << endl;

        cout << "     <list of input eventdisplay files (MC)> : test files with input evndisplay files" << endl;
        cout << "     <train vs test fraction> fraction of events to be used for training (typical 0.5)" << endl;
        cout << "     <reconstruction ID>:  e.g. 0,1,2,3" << endl;
        cout << "     telescope type ID (if not given: all telescope types are used)" << endl;
        cout << "                       (for VTS - these are telescope numbers)" << endl;
        cout << "     optional: train for energy/core reconstruction = \"BDTDispEnergy\"/\"BDTDispCore\"";
        cout << "(default = \"BDTDisp\": train for angular reconstrution)" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string       fInputFile  = argv[1] ;
    string       fOutputDir  = argv[2] ;
    float        fTrainTest  = atof( argv[3] );
    unsigned int iRecID      = atoi( argv[4] );
    ULong64_t    iTelType    = atoi( argv[5] ) ;
    string       iTargetBDT  = "BDTDisp";
    string iTMVAOptions = "VarTransform=N:NTrees=200:BoostType=AdaBoost:MaxDepth=8";
    string       iDataDirectory = "";
    string       iLayoutFile = "";
    string       iQualityCut = "size>1.&&ntubes>4.&&width>0.&&width<2.&&length>0.&&length<10.&&tgrad_x<100.*100.&&loss<0.20"; //mono cuts
    //iQualityCut = iQualityCut + "&&cross<20.0&&EHeight<100.&&Rcore<2000."; //stereo cuts
    if( argc >=  7 )
    {
        iTargetBDT =       argv[6]   ;
    }
    if( argc >=  8 )
    {
        iTMVAOptions = argv[7];
    }
    if( argc >= 9 )
    {
        iLayoutFile = argv[8];
    }
    if( argc >= 10 )
    {
        iDataDirectory = argv[9];
    }
    if( argc >= 11 )
    {
        iQualityCut = argv[10];
    }
    
    ///////////////////////////
    // print runparameters to screen
    cout << "trainTMVAforAngularReconstruction (" << VGlobalRunParameter::getEVNDISP_VERSION() << ")" << endl;
    cout << "------------------------------------" << endl;
    cout << endl;
    cout << "input file list with evendisplay containing the training events: " << endl;
    cout << fInputFile << endl;
    cout << endl;
    cout << "training/testing fraction: " << fTrainTest << endl;
    if( iTelType > 0 )
    {
        cout << "training for telescope type " << iTelType << endl;
    }
    else
    {
        cout << "training using data from all telescope types" << endl;
    }
    cout << endl;
    cout << "using events for reconstruction ID " << iRecID << endl;
    
    /////////////////////////
    if( fTrainTest <= 0.0 || fTrainTest >= 1.0 )
    {
        cout << endl;
        cout << "Error, 4th argument <train vs test fraction> = '" << fTrainTest << "' must fall in the range 0.0 < x < 1.0" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    
    ///////////////////////////
    // output file
    ostringstream iFileName;
    iFileName << fOutputDir << "/" << iTargetBDT;
    if( iTelType != 0 )
    {
        iFileName << "_" << iTelType;
    }
    iFileName << ".root";
    TFile iO( iFileName.str().c_str(), "RECREATE" );
    if( iO.IsZombie() )
    {
        cout << "Error creating output file: " << iFileName.str() << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    //////////////////////
    // fill training file
    if( iDataDirectory.size() == 0 && !writeTrainingFile( fInputFile, iTelType, iRecID, iLayoutFile ) )
    {
        cout << "error writing training file " << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else if( iDataDirectory.size() != 0 && !readTrainingFile( iTargetBDT, iTelType, iDataDirectory ) )
    {
        cout << "error reading training file " << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    //////////////////////
    // write training tree to output file
    iO.cd();
    map< ULong64_t, TTree* >::iterator fMapOfTrainingTree_iter;
    for( fMapOfTrainingTree_iter = fMapOfTrainingTree.begin(); fMapOfTrainingTree_iter != fMapOfTrainingTree.end(); ++fMapOfTrainingTree_iter )
    {
        if( fMapOfTrainingTree_iter->second )
        {
            cout << "\t writing training tree for telescope type " << fMapOfTrainingTree_iter->first;
            cout << " with " << fMapOfTrainingTree_iter->second->GetEntries() << " entries" << endl;
            fMapOfTrainingTree_iter->second->Write();
        }
    } 
    //////////////////////
    // train TMVA
    cout << "Number of telescope types: " << fMapOfTrainingTree.size() << endl;
    for( fMapOfTrainingTree_iter = fMapOfTrainingTree.begin(); fMapOfTrainingTree_iter != fMapOfTrainingTree.end(); ++fMapOfTrainingTree_iter )
    {
        trainTMVA( fOutputDir, fTrainTest, 
                fMapOfTrainingTree_iter->first, 
                fMapOfTrainingTree_iter->second, 
                iTargetBDT, iTMVAOptions, iQualityCut );
    }
    
    //////////////////////
    // close output file
    iO.Close();
    
}
