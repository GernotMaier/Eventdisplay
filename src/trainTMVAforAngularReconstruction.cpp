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
#include "TMVA/DataLoader.h"
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
#include "VSimpleStereoReconstructor.h"
#include "VUtilities.h"

using namespace std;

/////////////////////////////////////////////////////
// one tree per telescope type
map< ULong64_t, TTree* > fMapOfTrainingTree;
map< ULong64_t, unsigned int > fMapOfNTelescopeType;
/////////////////////////////////////////////////////

/*
    train TVMA method and write results into the corresponding directory

    one MVA per telescope type

    Allowed targets for training:

    BDTDisp
    BDTDispError
    BDTDispEnergy
    BDTDispCore

    Default is the BDT trainer (but other MLs like MLPs are allowed, but not well tested)

*/
bool trainTMVA( string iOutputDir, float iTrainTest,
                ULong64_t iTelType, TTree* iDataTree,
                string iTargetML, string iTMVAOptions,
                string iQualityCut, bool iSingleTelescopeAnalysis,
                bool iUseImageParameterErrors )
{
    cout << endl;
    cout << "Starting " << iTargetML;
    cout << " training for telescope type " << iTelType << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << endl;

    if(!iDataTree )
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
    ntrain = floor( nentries* iTrainTest ) ;
    ntest  = nentries - ntrain ;
    if( ntrain <= 100 || ntest <= 100 )
    {
        cout << endl;
        cout << "Error, train/test fraction is so small that only " << ntrain << "(" << ntest << ")";
        cout << " events were selected for training, while TMVA usually needs thousands of training/testing events to work properly.";
        cout << "Try increasing the train/test fraction... (you only have " << nentries;
        cout << " total events to designate for either training or testing...)" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }
    // unclear why factor of 0.8
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
    train_and_test_conditions << "V=True" << ":" ;
    train_and_test_conditions << "VerboseLevel=Info" << ":" ;
    train_and_test_conditions << "ScaleWithPreselEff=True";
    cout << "Train and test condition: " << train_and_test_conditions.str() << endl;
    cout << endl;

    // output file name
    ostringstream iFileName;
    iFileName << iOutputDir << "/" << iTargetML << "_" << iTelType << ".tmva.root";
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
    TMVA::Factory* factory = new TMVA::Factory( iTargetML.c_str(), i_tmva,
        "V:!DrawProgressBar:!Color:!Silent:AnalysisType=Regression:VerboseLevel=Debug:Correlations=True" );
    factory->SetVerbose( true );
    TMVA::DataLoader* dataloader = new TMVA::DataLoader( "" );

    // list of variables used by MVA method
    dataloader->AddVariable( "width", 'F' );
    dataloader->AddVariable( "length", 'F' );
    dataloader->AddVariable( "wol",    'F' );
    dataloader->AddVariable( "size", 'F' );
    dataloader->AddVariable( "ntubes", 'F' );
    // hard coded ASTRI telescope type
    // (no time gradient is available)
    if( iTelType != 201511619 )
    {
        dataloader->AddVariable( "tgrad_x*tgrad_x", 'F' );
    }
    if(!iSingleTelescopeAnalysis )
    {
        dataloader->AddVariable( "cross", 'F' );
    }
    dataloader->AddVariable( "asym", 'F' );
    dataloader->AddVariable( "loss", 'F' );
    dataloader->AddVariable( "dist", 'F' );
    dataloader->AddVariable( "fui", 'F' );
    if( iTargetML.find( "DispEnergy" ) != string::npos && !iSingleTelescopeAnalysis )
    {
        dataloader->AddVariable( "EHeight", 'F' );
        dataloader->AddVariable( "Rcore", 'F' );
    }
    if( iUseImageParameterErrors )
    {
        dataloader->AddVariable( "dispImageError", 'F' );
    }
    // spectators
    dataloader->AddSpectator( "cen_x", 'F' );
    dataloader->AddSpectator( "cen_y", 'F' );
    dataloader->AddSpectator( "cosphi", 'F' );
    dataloader->AddSpectator( "sinphi", 'F' );
    dataloader->AddSpectator( "MCe0", 'F' );
    dataloader->AddSpectator( "MCxoff", 'F' );
    dataloader->AddSpectator( "MCyoff", 'F' );
    dataloader->AddSpectator( "MCxcore", 'F' );
    dataloader->AddSpectator( "MCycore", 'F' );
    dataloader->AddSpectator( "MCrcore", 'F' );
    dataloader->AddSpectator( "NImages", 'F' );
    // train for energy reconstruction
    if( iTargetML.find( "DispEnergy" ) != string::npos )
    {
        dataloader->AddTarget( "dispEnergy", 'F' );
    }
    // rotation angle
    else if( iTargetML.find( "DispPhi" ) != string::npos )
    {
        dataloader->AddSpectator( "disp", 'F' );
        dataloader->AddSpectator( "dispError", 'F' );
        dataloader->AddTarget( "dispPhi", 'F' );
    }
    // train for error on disp reconstruction
    else if( iTargetML.find( "DispError" ) != string::npos )
    {
        dataloader->AddSpectator( "disp", 'F' );
        dataloader->AddSpectator( "dispPhi", 'F' );
        dataloader->AddTarget( "dispError", 'F', "dispError", 0., 10. );
    }
    // train for core reconstruction
    else if( iTargetML.find( "DispCore" ) != string::npos )
    {
        dataloader->AddTarget( "dispCore", 'F', "m", 0., 1.e5 );
    }
    // train for direction reconstruction
    else
    {
        dataloader->AddSpectator( "dispError", 'F' );
        dataloader->AddSpectator( "dispPhi", 'F' );
        dataloader->AddTarget( "disp", 'F' );
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
    if( iTargetML.find( "MLP" ) != string::npos )
    {
        iMVAName << "MLP_" << iTelType;
    }
    else
    {
        iMVAName << "BDT_" << iTelType;
    }
    sprintf( htitle, "%s", iMVAName.str().c_str() );

    TString methodstr( iTMVAOptions.c_str() );

    cout << "Built MethodStringStream: " << iTMVAOptions << endl;
    cout << endl;
    TString methodTitle( htitle );
    if( iTargetML.find( "MLP" ) != string::npos )
    {
        factory->BookMethod( dataloader, TMVA::Types::kMLP,  methodTitle, methodstr );
    }
    else
    {
        factory->BookMethod( dataloader, TMVA::Types::kBDT, methodTitle, methodstr ) ;
    }

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
    if(!is )
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
bool readTrainingFile( string iTargetML, ULong64_t iTelType, string iDataDirectory )
{
    fMapOfTrainingTree.clear();

    ostringstream iFileName;
    iFileName << iDataDirectory << "/" << iTargetML;
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
        if(!is )
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
                        unsigned int iRecID, string iArrayList,
                        bool redo_stereo_reconstruction )
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
    cout << "reading telescope list from ";
    cout << iInputFileList[0] << endl;

    Ctelconfig i_tel(&i_telChain );
    i_tel.GetEntry( 0 );
    unsigned int i_ntel = i_tel.NTel;

    vector< unsigned int > iHyperArrayID;
    vector< float > iFOV_tel;
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
    // (unfortunately inconsistent in data
    //  types required)
    vector< float > fTelX;
    vector< float > fTelY;
    vector< float > fTelZ;
    vector< double > fEM_TelX;
    vector< double > fEM_TelY;
    vector< double > fEM_TelZ;
    vector< ULong64_t > fTelType;
    unsigned int f_ntelType = 0;
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        i_tel.GetEntry( i );

        fTelX.push_back( i_tel.TelX );
        fTelY.push_back( i_tel.TelY );
        fTelZ.push_back( i_tel.TelZ );
        fEM_TelX.push_back( i_tel.TelX );
        fEM_TelY.push_back( i_tel.TelY );
        fEM_TelZ.push_back( i_tel.TelZ );
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
    float dcen_x = -1.;
    float dcen_y = -1.;
    float sinphi = -1.;
    float cosphi = -1.;
    float dphi = -1.;
    float size = -1.;    // actually log10(size)
    float ntubes = -1.;
    float loss = -1.;
    float asym = -1.;
    float width = -1.;
    float dwidth = -1.;
    float length = -1.;
    float dlength = -1.;
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
    float dispImageError = -1.;
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
    int   Fitstat = -1;

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

            fMapOfTrainingTree[i_tel.TelType]->Branch( "runNumber", &runNumber, "runNumber/I" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "eventNumber", &eventNumber, "eventNumber/I" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "tel", &tel, "tel/i" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cen_x", &cen_x, "cen_x/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cen_y", &cen_y, "cen_y/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "sinphi", &sinphi, "sinphi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cosphi", &cosphi, "cosphi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "size", &size, "size/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "ntubes", &ntubes, "ntubes/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "loss", &loss, "loss/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "asym", &asym, "asym/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "width", &width, "width/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "length", &length, "length/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "wol", &wol, "wol/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dist", &dist, "dist/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "fui", &fui, "fui/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "tgrad_x", &tgrad_x, "tgrad_x/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "meanPedvar_Image", &meanPedvar_Image, "meanPedvar_Image/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Fitstat", &Fitstat, "Fitstat/I" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dcen_x", &dcen_x, "dcen_x/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dcen_y", &dcen_y, "dcen_y/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dwidth", &dwidth, "dwidth/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dlength", &dlength, "dlength/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dphi", &dphi, "dphi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCe0", &MCe0, "MCe0/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCxoff", &MCxoff, "MCxoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCyoff", &MCyoff, "MCyoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCxcore", &MCxcore, "MCxcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCycore", &MCycore, "MCycore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCrcore", &MCrcore, "MCrcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Xcore", &Xcore, "Xcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Ycore", &Ycore, "Ycore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Rcore", &Rcore, "Rcore/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Xoff", &Xoff, "Xoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Yoff", &Yoff, "Yoff/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "LTrig", &LTrig, "LTrig/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "NImages", &NImages, "NImages/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "EHeight", &EmissionHeight, "EHeight/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCaz", &MCaz, "MCaz/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "MCze", &MCze, "MCze/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Ze", &ze, "Ze/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "Az", &az, "Az/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "disp", &disp, "disp/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispError", &dispError, "dispError/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispImageError", &dispImageError, "dispImageError/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "cross", &cross, "cross/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispPhi", &dispPhi, "dispPhi/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispEnergy", &dispEnergy, "dispEnergy/F" );
            fMapOfTrainingTree[i_tel.TelType]->Branch( "dispCore", &dispCore, "dispCore/F" );

        }
    }

    ////////////////////////////////////////////
    // filling of training trees;
    cout << "filling training trees for " << fMapOfTrainingTree.size() << " telescope type(s)" << endl;
    cout << "\t found " << f_ntelType << " telescopes of telescope type " << iTelType << endl;
    bool iSingleTelescopeAnalysis = false;
    if( f_ntelType == 1 )
    {
        iSingleTelescopeAnalysis = true;
        cout << "\t single telescope analysis" << endl;
    }
    fMapOfNTelescopeType[iTelType] = f_ntelType;

    // get showerpars tree
    TChain i_showerparsTree( "showerpars" );
    for( unsigned int f = 0; f < iInputFileList.size(); f++ )
    {
        i_showerparsTree.Add( iInputFileList[f].c_str(), 0 );
    }
    Cshowerpars i_showerpars(&i_showerparsTree, true, true );

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
    double fEM_width[fTelType.size()];
    double fEM_length[fTelType.size()];
    double fEM_cosphi[fTelType.size()];
    double fEM_sinphi[fTelType.size()];
    double fEM_weight[fTelType.size()];

    // stereo (intersection of line) reconstruction
    // needed for the re-calculation of 'cross'
    VSimpleStereoReconstructor i_SR;
    i_SR.initialize();

    /////////////////////////////////////////////////
    // loop over all events in trees
    int nentries = i_showerpars.fChain->GetEntries();
    cout << "Loop over " << nentries << " entries in source files" << endl;

    for( int n = 0; n < nentries; n++ )
    {
        // read events from event trees
        i_showerpars.GetEntry( n );

        // check recid
        if( iRecID >= i_showerpars.NMethods )
        {
            cout << "Error: invalid reconstruction ID.";
            cout << " Maximum allowed value is " << i_showerpars.NMethods << endl;
            return false;
        }

        // require:
        // - reconstructed event
        // - at least two telescopes
        if(!iSingleTelescopeAnalysis )
        {
            if( i_showerpars.Chi2[iRecID] < -999.
                    ||  i_showerpars.NImages[iRecID] < 2 )
            {
                continue;
            }
        }
        else
        {
            if( i_showerpars.NImages[iRecID] < 1 )
            {
                continue;
            }
        }

        // check if there are image of this telescope type
        // (hyper-array)
        int i_nteltypecounter = 0;
        for( unsigned int i = 0; i < fTelType.size(); i++ )
        {
            if(( fTelType[i] == iTelType || iTelType == 0 )
                    && ( int )i_showerpars.ImgSel_list[iRecID][i] > 0 )
            {
                i_nteltypecounter++;
            }
        }
        if( i_nteltypecounter == 0 )
        {
            continue;
        }

        /////////////////////////////////////////////////////////
        // calculate emission height and cross
        for( unsigned int i = 0; i < i_tpars.size(); i++ )
        {
            fEM_size[i] = -1.;
            fEM_cen_x[i] = 0.;
            fEM_cen_y[i] = 0.;
            fEM_width[i] = 0.;
            fEM_length[i] = 0.;
            fEM_cosphi[i] = 0.;
            fEM_sinphi[i] = 0.;
            fEM_weight[i] = 0.;

            if(( int )i_showerpars.ImgSel_list[iRecID][i] < 1
                    && ( i_showerpars.NImages[iRecID] > 1 || !iSingleTelescopeAnalysis ) )
            {
                continue;
            }
            // note: this is different to what happens in the analysis:
            // here the emissionheight / direction is calculated only from the
            // telescopes of the given telescope type
            if(!i_tpars[i] )
            {
                continue;
            }

            i_tpars[i]->GetEntry( n );

            if( i_tpars[i]->size > 0. )
            {
                fEM_size[i] = i_tpars[i]->size;
                fEM_cen_x[i] = i_tpars[i]->cen_x;
                fEM_cen_y[i] = i_tpars[i]->cen_y;
                fEM_width[i] = i_tpars[i]->width;
                fEM_length[i] = i_tpars[i]->length;
                fEM_cosphi[i] = i_tpars[i]->cosphi;
                fEM_sinphi[i] = i_tpars[i]->sinphi;
                // weight is always 1: all telescopes
                // are of same type
                fEM_weight[i] = 1.;
            }
        }
        EmissionHeight = fEmissionHeightCalculator->getEmissionHeight( fEM_cen_x, fEM_cen_y, fEM_size,
                         i_showerpars.ArrayPointing_Azimuth,
                         i_showerpars.ArrayPointing_Elevation );

        if(!iSingleTelescopeAnalysis )
        {
            i_SR.reconstruct_direction_and_core(
                fEM_TelX.size(),
                i_showerpars.ArrayPointing_Elevation,
                i_showerpars.ArrayPointing_Azimuth,
                &fEM_TelX[0], &fEM_TelY[0], &fEM_TelZ[0],
                fEM_size,
                fEM_cen_x,
                fEM_cen_y,
                fEM_cosphi,
                fEM_sinphi,
                fEM_width,
                fEM_length,
                fEM_weight );
        }

        //////////////////////////////////////
        // loop over all telescopes
        for( unsigned int i = 0; i < i_tpars.size(); i++ )
        {
            // check if telescope was reconstructed
            if(( int )i_showerpars.ImgSel_list[iRecID][i] < 1 )
            {
                continue;
            }
            // check if telescope is of valid telescope type
            if(( fTelType[i] != iTelType && iTelType != 0 ) || !i_tpars[i] )
            {
                continue;
            }
            // check if event is not completely out of the FOV
            // (use 20% x size of the camera)
            if( i < iFOV_tel.size()
                    && sqrt( i_showerpars.MCxoff * i_showerpars.MCxoff
                        + i_showerpars.MCyoff * i_showerpars.MCyoff ) > iFOV_tel[i] * 0.5 * 1.2 )
            {
                continue;
            }
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
            ntubes      = log10( i_tpars[i]->ntubes );
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
            Fitstat     = i_tpars[i]->Fitstat;
            if( i_tpars[i]->hasParameterErrors() )
            {
                dcen_x = i_tpars[i]->dcen_x;
                dcen_y = i_tpars[i]->dcen_y;
                dwidth = i_tpars[i]->dwidth;
                dlength = i_tpars[i]->dlength;
                dphi = i_tpars[i]->dphi;
            }
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

            Rcore       = VUtilities::line_point_distance( Ycore, -1.*Xcore,   0., ze, az, fTelY[i], -1.*fTelX[i], fTelZ[i] );
            MCrcore     = VUtilities::line_point_distance( MCycore, -1.*MCxcore, 0., MCze, MCaz, fTelY[i], -1.*fTelX[i], fTelZ[i] );

            if( Rcore < 0. && !iSingleTelescopeAnalysis )
            {
                continue;
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////
            // calculate disp (observe sign convention for MC in y direction for MCyoff and Yoff)
            disp  = sqrt(( cen_y + MCyoff ) * ( cen_y + MCyoff ) + ( cen_x - MCxoff ) * ( cen_x - MCxoff ) );
            if( iSingleTelescopeAnalysis )
            {
                cross = 0.;
            }
            else if( redo_stereo_reconstruction )
            {
                cross = sqrt(( cen_y + i_SR.fShower_Yoffset ) * ( cen_y + i_SR.fShower_Yoffset )
                             + ( cen_x - i_SR.fShower_Xoffset ) * ( cen_x - i_SR.fShower_Xoffset ) );
            }
            else
            {
                cross = sqrt(( cen_y + Yoff ) * ( cen_y + Yoff )
                             + ( cen_x - Xoff ) * ( cen_x - Xoff ) );
            }
            dispPhi = TMath::ATan2( sinphi, cosphi ) - TMath::ATan2( cen_y + MCyoff, cen_x - MCxoff );

            // disp error: the expected difference between true and
            //             reconstructed direction
            // Note that this is only the error expected due to the mismatch
            // of the image length axis with the true direction, not the error
            // due to a wrong prediction of disp by the BDT
            dispError = 0;
            float x1 = cen_x - disp * cosphi;
            float x2 = cen_x + disp * cosphi;
            float y1 = cen_y - disp * sinphi;
            float y2 = cen_y + disp * sinphi;
            if( sqrt(( x1 - MCxoff ) * ( x1 - MCxoff ) + ( y1 + MCyoff ) * ( y1 + MCyoff ) )
                    < sqrt(( x2 - MCxoff ) * ( x2 - MCxoff ) + ( y2 + MCyoff ) * ( y2 + MCyoff ) ) )
            {
                dispError = sqrt(( x1 - MCxoff ) * ( x1 - MCxoff ) + ( y1 + MCyoff ) * ( y1 + MCyoff ) );
            }
            else
            {
                dispError = sqrt(( x2 - MCxoff ) * ( x2 - MCxoff ) + ( y2 + MCyoff ) * ( y2 + MCyoff ) );
            }
            // disp uncertainty
            // - only possible for LL image fitting
            // - standard error propagation for disp calculation
            // - assume LL fit is performed in rotated frame
            //   (meaning: cen_x and length on one axis)
            dispImageError = -1.;
            if( i_tpars[i]->hasParameterErrors() )
            {
                dispImageError = sqrt( dcen_x* dcen_x
                                       +  dlength* dlength );
            }

            dispEnergy = i_showerpars.MCe0;
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
        cout << "./trainTMVAforAngularReconstruction <list of input eventdisplay files (MC)> <output directory>" << endl;
        cout << "                                     <train vs test fraction> <RecID> <telescope type>" << endl;
        cout << "                                     [train for angular / energy / core reconstruction]" << endl;
        cout << "                                     [MVA options] [array layout file] [directory with training trees]" << endl;
        cout << "                                     [quality cut] [use image parameter errors (default=off=0)]";
        cout << endl;
        cout << endl;

        cout << "     <list of input eventdisplay files (MC)> : test files with input evndisplay files" << endl;
        cout << "     <train vs test fraction> fraction of events to be used for training (typical 0.5)" << endl;
        cout << "     <reconstruction ID>:  e.g. 0,1,2,3" << endl;
        cout << "     telescope type ID (if not given: all telescope types are used)" << endl;
        cout << "                       (for VTS - these are telescope numbers)" << endl;
        cout << "     optional: train for energy/core reconstruction = \"BDTDispEnergy\"/\"BDTDispCore\"";
        cout << "(default = \"BDTDisp\": train for angular reconstruction)" << endl;
        cout << endl;
        exit( EXIT_SUCCESS );
    }
    string       fInputFile  = argv[1] ;
    string       fOutputDir  = argv[2] ;
    float        fTrainTest  = atof( argv[3] );
    unsigned int iRecID      = atoi( argv[4] );
    ULong64_t    iTelType    = atoi( argv[5] ) ;
    string       iTargetML  = "BDTDisp";
    if( argc >=  7 )
    {
        iTargetML = argv[6];
    }
    // tmva options likely overwritten from command line
    string iTMVAOptions = "VarTransform=N:NTrees=200:BoostType=AdaBoost:MaxDepth=8";
    if( argc >=  8 )
    {
        iTMVAOptions = argv[7];
    }
    string       iLayoutFile = "";
    if( argc >= 9 )
    {
        iLayoutFile = argv[8];
    }
    string       iDataDirectory = "";
    if( argc >= 10 )
    {
        iDataDirectory = argv[9];
    }
    // quality cut likely overwritten from command line
    string       iQualityCut = "size>1.&&ntubes>log10(4.)&&width>0.&&width<2.&&length>0.&&length<10.";
    iQualityCut = iQualityCut + "&&tgrad_x<100.*100.&&loss<0.20&&cross<20.0&&Rcore<2000.";
    if( argc >= 11 )
    {
        iQualityCut = argv[10];
    }
    bool iUseImageParameterErrors = false;
    if( argc >= 12 )
    {
        iUseImageParameterErrors = ( bool )( atoi( argv[11] ) );
    }
    bool redo_stereo_reconstruction = false;

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
    if( iUseImageParameterErrors )
    {
        cout << "using image parameter fit errors" << endl;
    }

    /////////////////////////
    if( fTrainTest <= 0.0 || fTrainTest >= 1.0 )
    {
        cout << endl;
        cout << "Error, <train vs test fraction> = '";
        cout << fTrainTest << "' must fall in the range 0.0 < x < 1.0" << endl;
        cout << endl;
        exit( EXIT_FAILURE );
    }

    ///////////////////////////
    // output file
    ostringstream iFileName;
    iFileName << fOutputDir << "/" << iTargetML;
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
    if( iDataDirectory.size() == 0
            && !writeTrainingFile( fInputFile, iTelType, iRecID, iLayoutFile, redo_stereo_reconstruction ) )
    {
        cout << "error writing training file " << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    else if( iDataDirectory.size() != 0
            && !readTrainingFile( iTargetML, iTelType, iDataDirectory ) )
    {
        cout << "error reading training file " << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    //////////////////////
    // write training tree to output file
    iO.cd();
    map< ULong64_t, TTree* >::iterator fMapOfTrainingTree_iter;
    for( fMapOfTrainingTree_iter = fMapOfTrainingTree.begin();
            fMapOfTrainingTree_iter != fMapOfTrainingTree.end();
            ++fMapOfTrainingTree_iter )
    {
        if( fMapOfTrainingTree_iter->second )
        {
            cout << "\t writing training tree for telescope type " << fMapOfTrainingTree_iter->first;
            cout << " with " << fMapOfTrainingTree_iter->second->GetEntries() << " entries " << endl;
            cout << "to " << iFileName.str() << endl;
            fMapOfTrainingTree_iter->second->Write();
        }
    }
    //////////////////////
    // train TMVA
    cout << "Number of telescope types: " << fMapOfTrainingTree.size() << endl;
    for( fMapOfTrainingTree_iter = fMapOfTrainingTree.begin();
            fMapOfTrainingTree_iter != fMapOfTrainingTree.end();
            ++fMapOfTrainingTree_iter )
    {
        bool iSingleTel = false;
        if( fMapOfNTelescopeType.find( fMapOfTrainingTree_iter->first ) != fMapOfNTelescopeType.end() )
        {
            if( fMapOfNTelescopeType[fMapOfTrainingTree_iter->first] == 1 )
            {
                iSingleTel = true;
            }
        }
        trainTMVA( fOutputDir, fTrainTest,
                   fMapOfTrainingTree_iter->first,
                   fMapOfTrainingTree_iter->second,
                   iTargetML, iTMVAOptions, iQualityCut,
                   iSingleTel, iUseImageParameterErrors );
    }

    //////////////////////
    // close output file
    iO.Close();

}
