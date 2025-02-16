//! VTableLookupDataHandler data class for mscw and energy reconstruction

#ifndef VTableLookupDataHandler_H
#define VTableLookupDataHandler_H

#include "VEmissionHeightCalculator.h"
#include "VDeadTime.h"
#include "VEffectiveAreaCalculatorMCHistograms.h"
#include "VMonteCarloRunHeader.h"
#include "VDispAnalyzer.h"
#include "VPointingCorrectionsTreeReader.h"
#include "VSimpleStereoReconstructor.h"
#include "VTableLookupRunParameter.h"
#include "VUtilities.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TKey.h"

#include "Cshowerpars.h"
#include "Ctelconfig.h"
#include "Ctpars.h"

#include <bitset>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class VTableLookupDataHandler
{
    private:
        unsigned int fDebug;

        VTableLookupRunParameter* fTLRunParameter;  //! lookup parameters

        UInt_t fNMethods;                            //!< number of direction and core reconstruction methods in evndisp
        vector< string > finputfile;              //!< input file name
        string foutputfile;                       //!< output file name
        TFile* fOutFile;                          //!< point to output file
        bool   fwrite;                            //!< true for table filling

        unsigned int fNTel;                       //!< number of telescopes
        unsigned int fNTelComb;                   //!< number of telescope combinations
        Long64_t fNEntries;                       //!< total number of events in input tree
        int fEventCounter;                        //!< event counter for input tree
        int fMethod;                              //!< which direction and core reconstruction method data should be used

        double fMaxTotalTime;                     //!< time to analyse in a run (in [s])
        double fTotalTime;                        //!< time passed in analysing the run (in [s])
        double fTotalTime0;                       //!< time of first event (in [s])

        VEmissionHeightCalculator* fEmissionHeightCalculator;
        VDispAnalyzer*             fDispAnalyzerDirection;
        VDispAnalyzer*             fDispAnalyzerEnergy;
        VDispAnalyzer*             fDispAnalyzerCore;
        VDispAnalyzer*             fDispAnalyzerDirectionError;

        double fSelectRandom;
        int fSelectRandomSeed;
        TRandom3* fRandom;

        int fSSR_NImages_min;
        float fSSR_AxesAngles_min;

        // MC parameter
        unsigned int fMinImages;
        double fMC_distance_to_cameracenter_min;
        double fMC_distance_to_cameracenter_max;

        double fMCSpectralIndex;
        double fMCMinEnergy;
        double fMCMaxEnergy;
        double fSpectralIndex;

        vector< double > fNoiseLevel;
        vector< double > fCurrentNoiseLevel;

        // input trees
        int fEventDisplayFileFormat;
        TChain* fTshowerpars;
        TChain* fTshowerpars_QCCut;
        Cshowerpars* fshowerpars;
        TChain* fTtelconfig;
        Ctelconfig* ftelconfig;
        vector< TChain* > fTtpars;
        vector< Ctpars* > ftpars;
        vector< VPointingCorrectionsTreeReader* > fpointingCorrections;
        TChain* fDeepLearnerpars;

        double fEventWeight;

        // MC energy histograms
        TH1D* hE0mc;
        TH2D* hXYmc;
        TH2D* hDE0mc;
        TH2D* hWE0mc;
        TH1D* hZe;
        TH1D* hE0trig;
        TH2D* hXYtrig;
        TH2D* hDE0trig;
        TH2D* hWE0trig;
        TList* hisList;

        // telescope positions from fTtelconfig
        double fTelX[VDST_MAXTELESCOPES];
        double fTelY[VDST_MAXTELESCOPES];
        double fTelZ[VDST_MAXTELESCOPES];
        double fFocalLength[VDST_MAXTELESCOPES];
        double fTelFOV[VDST_MAXTELESCOPES];
        ULong64_t fTel_type[VDST_MAXTELESCOPES];
        map<ULong64_t, unsigned int > fList_of_Tel_type;                      // [teltype][number of telescopes for this type]
        map<ULong64_t, unsigned int >::iterator fList_of_Tel_type_iterator;
        vector< unsigned int > fTel_type_counter;
        // telescope pointing
        float  fArrayElevation;
        float  fArrayAzimuth;
        double fTelElevation[VDST_MAXTELESCOPES];
        double fTelAzimuth[VDST_MAXTELESCOPES];
        double fWobbleN;
        double fWobbleE;
        float  fArrayPointing_Elevation;
        float  fArrayPointing_Azimuth;
        float  fArrayPointing_RotationAngle;

        // output trees
        TTree* fOTree;
        bool fShortTree;                          //!< use short version of output tree
        bool bWriteMCPars;
        bool fTreeWithParameterErrors;

        // cut statistics
        bool fEventStatus;
        unsigned int fNStats_All;
        unsigned int fNStats_Rec;
        unsigned int fNStats_NImagesCut;
        unsigned int fNStats_Chi2Cut;
        unsigned int fNStats_CoreErrorCut;
        unsigned int fNStats_WobbleCut;
        unsigned int fNStats_WobbleMinCut;
        unsigned int fNStats_WobbleMaxCut;

        void   calcDistances();                //!< calculate distances between telescopes and shower core
        void   calcEmissionHeights();
        double calculateMeanNoiseLevel( bool bCurrentNoiseLevel = false );
        bool   checkIfFilesInChainAreRecovered( TChain* c );
        void   copyMCHistograms();
        bool   copyMCRunheader();
        void   copyTree_from_evndispFile( string iTreename = "MCpars" );
        void   copy_telconfig();
        bool   doImageQualitySelection( unsigned int iTelID );
        void   doStereoReconstruction();
        void   initializeTelTypeVector();
        int    fillNextEvent( bool bShort );
        void   printCutStatistics();
        void   printTelescopesList( unsigned int iPrintParameter );
        bool   randomSelected();
        void   resetImageParameters();
        void   resetImageParameters( unsigned int i );
        void   setEventWeightfromMCSpectrum();
        void   setSelectRandom( double iX, int iS );
        void   writeDeadTimeHistograms();

    public:

        //  data written to output file

        int runNumber;
        int eventNumber;
        int MJD;
        double time;

        bool fIsMC;                               //!< data is MC
        int    fMCPrimary;
        double fMCEnergy;                         //!< MC energy
        double fMCxcore;
        double fMCycore;
        double fMCxcore_SC;
        double fMCycore_SC;
        double fMCxcos;
        double fMCycos;
        double fMCaz;
        double fMCze;
        double fMCxoff;
        double fMCyoff;
        Int_t	       fMCCorsikaRunID;
        Int_t	       fMCCorsikaShowerID;
        Float_t        fMCFirstInteractionHeight;
        Float_t	       fMCFirstInteractionDepth;

        ULong64_t LTrig;
        unsigned int fNTrig;
        int fNImages;
        UShort_t fNImages_QCTree[VDST_MAXTELESCOPES];
        ULong64_t fImgSel;
        bool fImgSel_list[VDST_MAXTELESCOPES];                   // length of array: number of telescopes
        unsigned int fImgSel_list_short[VDST_MAXTELESCOPES];     // length of array: number of active telescopes
        int          fNTelTypes;
        ULong64_t    fTtypeID[VDST_MAXTELESCOPES];                // [tel type counter] = teltype
        unsigned int NImages_Ttype[VDST_MAXTELESCOPES];
        double fimg2_ang;
        double fZe;                               //!< zenith angle
        double fAz;
        double fRA;
        double fDec;
        double fXoff;
        double fYoff;
        double fXoff_derot;
        double fYoff_derot;
        double fstdS;
        double ftheta2;
        double fXcore;
        double fYcore;
        double fXcore_SC;
        double fYcore_SC;
        double fstdP;
        double fchi2;                             //!< chi2 from array reconstruction
        Float_t fchi2_QCTree[VDST_MAXTELESCOPES];
        //
        // {1}
        double fMCEnergyArray [VDST_MAXTELESCOPES];
        float  fmeanPedvar_ImageT[VDST_MAXTELESCOPES];
        float  fmeanPedvar_Image;
        double fdist     [VDST_MAXTELESCOPES];
        double ffui       [VDST_MAXTELESCOPES];
        double fdist_telType[VDST_MAXTELESCOPES];
        double fsize     [VDST_MAXTELESCOPES];
        double fsize2    [VDST_MAXTELESCOPES];
        double fsizeCorr [VDST_MAXTELESCOPES];
        double fsize_telType[VDST_MAXTELESCOPES];
        double floss     [VDST_MAXTELESCOPES];
        double fLoss_telType[VDST_MAXTELESCOPES];
        double fDistance_telType[VDST_MAXTELESCOPES];
        double ffracLow  [VDST_MAXTELESCOPES];
        double fmax1     [VDST_MAXTELESCOPES];
        double fmax2     [VDST_MAXTELESCOPES];
        double fmax3     [VDST_MAXTELESCOPES];
        int    fmaxindex1     [VDST_MAXTELESCOPES];
        int    fmaxindex2     [VDST_MAXTELESCOPES];
        int    fmaxindex3     [VDST_MAXTELESCOPES];
        double fwidth    [VDST_MAXTELESCOPES];
        double fdwidth    [VDST_MAXTELESCOPES];
        double fwidth_telType[VDST_MAXTELESCOPES];
        double flength   [VDST_MAXTELESCOPES];
        double fdlength   [VDST_MAXTELESCOPES];
        double flength_telType[VDST_MAXTELESCOPES];
        int    fntubes   [VDST_MAXTELESCOPES];
        unsigned short int fnsat[VDST_MAXTELESCOPES];
        unsigned short int fnlowgain[VDST_MAXTELESCOPES];
        double falpha    [VDST_MAXTELESCOPES];
        double flos      [VDST_MAXTELESCOPES];
        double fasym     [VDST_MAXTELESCOPES];
        double fcen_x    [VDST_MAXTELESCOPES];
        double fdcen_x    [VDST_MAXTELESCOPES];
        double fcen_y    [VDST_MAXTELESCOPES];
        double fdcen_y    [VDST_MAXTELESCOPES];
        double fdphi[VDST_MAXTELESCOPES];
        double fcosphi   [VDST_MAXTELESCOPES];
        double fsinphi   [VDST_MAXTELESCOPES];
        double ftgrad_x  [VDST_MAXTELESCOPES];
        double ftgrad_x_telType[VDST_MAXTELESCOPES];
        double ftchisq_x [VDST_MAXTELESCOPES];
        double fweight   [VDST_MAXTELESCOPES];    //!< always 1.
        double fweightDispBDTs[VDST_MAXTELESCOPES];
        int    fFitstat  [VDST_MAXTELESCOPES];
        double fR        [VDST_MAXTELESCOPES];    //!< distance from each telescope to reconstructed shower core
        double fRTel        [VDST_MAXTELESCOPES];    //!< distance from each telescope to reconstructed shower core
        double fR_telType[VDST_MAXTELESCOPES];    //!< distance from each telescope to reconstructed shower core (depending on tel type)
        double fR_short[VDST_MAXTELESCOPES];    //!< distance from each telescope to reconstructed shower core (depending on imgsel)
        double fR_MC      [VDST_MAXTELESCOPES];    //!< distance from each telescope to MC shower core
        double fR_telType_MC[VDST_MAXTELESCOPES];    //!< distance from each telescope to MC shower core (depending on tel type)
        double fR_short_MC[VDST_MAXTELESCOPES];    //!< distance from each telescope to MC shower core (depending on imgsel)
        int    fntubes_short[VDST_MAXTELESCOPES];
        float  fdist_short[VDST_MAXTELESCOPES];
        float  fwidth_short[VDST_MAXTELESCOPES];
        float  fdwidth_short[VDST_MAXTELESCOPES];
        float  fasym_short[VDST_MAXTELESCOPES];
        float  flength_short[VDST_MAXTELESCOPES];
        float  fdlength_short[VDST_MAXTELESCOPES];
        float  fcen_x_short[VDST_MAXTELESCOPES];
        float  fcen_y_short[VDST_MAXTELESCOPES];
        float  fdcen_x_short[VDST_MAXTELESCOPES];
        float  fdcen_y_short[VDST_MAXTELESCOPES];
        float  fcosphi_short[VDST_MAXTELESCOPES];
        float  fsinphi_short[VDST_MAXTELESCOPES];
        float  fdphi_short[VDST_MAXTELESCOPES];
        float  floss_short[VDST_MAXTELESCOPES];
        float  ffui_short[VDST_MAXTELESCOPES];
        float  fsize_short[VDST_MAXTELESCOPES];
        float  fcross_short[VDST_MAXTELESCOPES];
        float  fcrossO_short[VDST_MAXTELESCOPES];
        float  ftgrad_x_short[VDST_MAXTELESCOPES];
        int    fFitstat_short[VDST_MAXTELESCOPES];

        UInt_t PixelListN[VDST_MAXTELESCOPES];
        UInt_t PixelListNPixelNN;
        UInt_t PixelID[VDST_MAXTELESCOPES* VDST_MAXCHANNELS];
        UInt_t PixelType[VDST_MAXTELESCOPES* VDST_MAXCHANNELS];
        Float_t PixelIntensity[VDST_MAXTELESCOPES* VDST_MAXCHANNELS];
        Float_t PixelTimingT0[VDST_MAXTELESCOPES* VDST_MAXCHANNELS];
        Float_t PixelPE[VDST_MAXTELESCOPES* VDST_MAXCHANNELS];

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        // results
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////

        int    fnmscw;                            //!< number of images used for mscw/mscl calculation
        int    fnenergyT;                         //!< number of images used for the energy calculation
        int    fenergyQL;                         //!< quality label for energy calculation
        // MSCW
        double fmscw;                             //!< mean scaled width
        double ftmscw    [VDST_MAXTELESCOPES];    //!< mscw assigned to each telescope
        float  ftmscw_sigma[VDST_MAXTELESCOPES];  //!< mscw  sigma  assigned to each telescope
        // MWR
        float  fmwr;                              //!< mean width ratio

        // MSCL
        double fmscl;                             //!< mean scaled length
        double ftmscl    [VDST_MAXTELESCOPES];    //!< mscl assigned to each telescope
        float  ftmscl_sigma[VDST_MAXTELESCOPES];  //!< mscl  sigma  assigned to each telescope
        // MWL
        float  fmlr;                              //!< mean length ratio

        // reconstructed energy
        double fenergyS;                          //!< reconstructed primary energy
        double fechi2S;                           //!< chi2 from reconstructed primary energy
        double fdES;                              //!< dE from reconstructed primary energy
        double feAbsError;                        //!< median absolute error of reconstructed energy
        double fES       [VDST_MAXTELESCOPES];    //!< energy assigned to each telescope (method 1)
        double fES_short [VDST_MAXTELESCOPES];    //!< energy assigned to each telescope (method 1) (depending on imgsel)

        // time gradient
        double fmsct;                             //!< mean scaled time gradient
        double ftmsct    [VDST_MAXTELESCOPES];    //!< msct assigned to each telescope
        float  ftmsct_sigma[VDST_MAXTELESCOPES];  //!< msct sigma assigned to each telescope

        // emission height
        unsigned int fNTelPairs;
        float  fEmissionHeightMean;
        float  fEmissionHeightChi2;
        //!< (note that maximum array length should be larger than MaxNbrTel
        float  fEmissionHeightT[VDST_MAXTELESCOPES];

        double fSizeSecondMax;
        double ftheta2_All[25];

        // disp related variables
        float fXoff_edisp;
        float fYoff_edisp;
        float fXoff_intersect;                  //! keep direction from intersection method
        float fYoff_intersect;                  //! keep direction from intersection method
        float fXoff_T[VDST_MAXTELESCOPES];      //! direction reconstructed for each telescope
        float fYoff_T[VDST_MAXTELESCOPES];      //! direction reconstructed for each telescope
        float fWoff_T[VDST_MAXTELESCOPES];      //! direction reconstructed for each telescope (weight)
        float fDoff_T[VDST_MAXTELESCOPES];      //! (disp value)
        float fDispAbsSumWeigth;                //! sum of absolute values of disp weights
        unsigned int fToff_T[VDST_MAXTELESCOPES]; //! list of telescope participating in disp
        unsigned int fnxyoff;                   //! number of images used for disp direction reconstruction
        // difference in disp event direction between telescopes
        double fDispDiff;

        // dispBDT cut parameter
        double fmaxdist_qc[VDST_MAXTELESCOPES];

        // deep learner parameters
        bool fIsDeepLearner;
        double dl_gammaness;
        Bool_t dl_isGamma;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        VTableLookupDataHandler( bool iWrite, VTableLookupRunParameter* iT = 0 );
        ~VTableLookupDataHandler() {}

        bool cut()                                //!< apply cuts on successful reconstruction to input data
        {
            return cut( false );
        }
        bool cut( bool bWrite );
        void fill();                              //!< fill output tree
        void fillMChistograms();
        void setfillTables( bool ib )
        {
            fwrite = ib;
        }
        int getArrayReconstructionMethod()
        {
            return fMethod;
        }
        double* getDistance()
        {
            return fdist;
        }
        double* getDistance( ULong64_t iTelType );
        double* getDistanceToCore()
        {
            return fR;
        }
        double* getDistanceToCoreTel()
        {
            return fRTel;
        }
        double* getDistanceToCore( ULong64_t iTelType, bool iMC = false );
        int    getEventNumber()
        {
            return eventNumber;
        }
        bool   getEventStatus()
        {
            return fEventStatus;
        }
        double getEventWeight()
        {
            return fEventWeight;
        }
        double getMCAz()
        {
            return fMCaz;
        }
        double getMCZe()
        {
            return fMCze;
        }
        double getAz()
        {
            return fAz;
        };
        double getMCDistance();
        double* getMCEnergyArray();
        double getMCEnergy()
        {
            return fMCEnergy;
        }
        double getMCWobbleOffset()
        {
            return sqrt( fMCxoff * fMCxoff + fMCyoff * fMCyoff );
        }
        int    getEventCounter()
        {
            return fEventCounter;
        }
        double* getLength()
        {
            return flength;
        }
        double* getLength( ULong64_t iTelType );
        map<ULong64_t, unsigned int > getList_of_Tel_type()
        {
            return fList_of_Tel_type;
        }
        double* getLoss()
        {
            return floss;
        }
        double* getLoss( ULong64_t iTelType );
        unsigned int getNTel_type( ULong64_t t )
        {
            if( fList_of_Tel_type.find( t ) != fList_of_Tel_type.end() )
            {
                return fList_of_Tel_type[t];
            }
            else
            {
                return 99999;
            }
        }
        unsigned int getMaxNbrTel() const
        {
            return VDST_MAXTELESCOPES;
        }
        int* getNtubes()
        {
            return fntubes;
        }
        double* getMSCWtel()
        {
            return ftmscw;
        }
        double* getMSCLtel()
        {
            return ftmscl;
        }
        double getMCMinEnergy()
        {
            return fMCMinEnergy;
        }
        double getMCMaxEnergy()
        {
            return fMCMaxEnergy;
        }
        double getMCSpectralIndex()
        {
            return fMCSpectralIndex;
        }
        double getSpectralIndex()
        {
            return fSpectralIndex;
        }
        Long64_t getNEntries()
        {
            return fNEntries;
        }
        int  getNEvents()
        {
            if( fTshowerpars )
            {
                return fTshowerpars->GetEntries();
            }
            else
            {
                return 0;
            }
        }
        bool getNextEvent( bool bShort = false );         //!< get next event from evndisp tree
        double getMeanNoiseLevel( bool bCurrentNoiseLevel = false )
        {
            if( !bCurrentNoiseLevel )
            {
                return calculateMeanNoiseLevel( false );
            }
            else
            {
                return ( double )fmeanPedvar_Image;
            }
        }
        vector< double > getNoiseLevel( bool bCurrentNoiseLevel = false )
        {
            if( !bCurrentNoiseLevel )
            {
                return fNoiseLevel;
            }
            else
            {
                return fCurrentNoiseLevel;
            }
        }
        map< ULong64_t, double > getNoiseLevel_per_TelescopeType();
        unsigned int getNTel()
        {
            return fNTel;
        }
        unsigned int getNTelTypes()
        {
            return fNTelTypes;
        }
        TFile* getOutputFile()
        {
            return fOutFile;
        }
        double* getSize( double iSizeCorrection = 1., bool iSelectedImagesOnly = false, bool iSize2 = false );
        double* getSize( double iSizeCorrection, ULong64_t iTelType, bool iSelectedImagesOnly, bool iSize2 = false );
        double* getSize2( double iSizeCorrection = 1., bool iSelectedImagesOnly = false )
        {
            return getSize( iSizeCorrection, iSelectedImagesOnly, true );
        }
        double* getSize2( double iSizeCorrection, ULong64_t iTelType, bool iSelectedImagesOnly )
        {
            return getSize( iSizeCorrection, iTelType, iSelectedImagesOnly, true );
        }
        double* getTimeGradient()
        {
            return ftgrad_x;
        }
        double* getTimeGradient( ULong64_t iTelType );
        double* getWeight(bool dispBDT=false)
        {
            if( dispBDT )
            {
                return fweightDispBDTs;
            }
            else
            {
                return fweight;
            }
        }
        double* getWidth()
        {
            return fwidth;
        }
        double* getWidth( ULong64_t iTelType );
        double getMaxTotalTime()
        {
            return fMaxTotalTime;
        }
        double getTelElevation( unsigned int iTelID )
        {
            return fTelElevation[iTelID];
        }
        double getTelElevation();
        ULong64_t getTelType( unsigned int iTelID )
        {
            return fTel_type[iTelID];
        }
        unsigned int getTelType_arraycounter( unsigned int iTelID );      // return position of tel type for this telescope in array counter
        double getZe();
        double getTheta2()
        {
            return ftheta2;
        }
        double getWobbleOffset()
        {
            if( fXoff < -9998. || fYoff < -9998. )
            {
                return -1.;
            }
            return sqrt( fXoff * fXoff + fYoff * fYoff );
        }
        bool isReconstructed( bool isReconstructed = false );
        bool readRunParameter();
        void reset();                             //!< reset a few output variables
        void resetAll();
        void setEnergy( double iES, double iChi2S, double idES, double ieAbsError = -999. )
        {
            fenergyS = iES;
            fechi2S = iChi2S;
            fdES = idES;
            feAbsError = ieAbsError;
        }
        void setEnergyT( int i, double iETS )
        {
            fES[i] = iETS;
        }
        bool setInputFile( vector< string > );              //!< set input file
        void setMCMinEnergy( double iB )
        {
            fMCMinEnergy = iB;
        }
        void setMCMaxEnergy( double iB )
        {
            fMCMaxEnergy = iB;
        }
        void setMCSpectralIndex( double iB )
        {
            fMCSpectralIndex = iB;
        }
        void setSpectralIndex( double iB )
        {
            fSpectralIndex = iB;
        }
        void setMethod( int );                    //!< data from which evndisp reconstruction method should be used?
        void setMSCL( double iMSLC )
        {
            fmscl = iMSLC;
        }
        void setMSCLT( int i, double iMSLC, float iMSLC_T = -99. )
        {
            ftmscl[i] = iMSLC;
            ftmscl_sigma[i] = iMSLC_T;
        }
        void setNEnergyT( int in )
        {
            fnenergyT = in;
        }
        void setNEnergyQuality( int in )
        {
            fenergyQL = in;
        }
        void setNMSCW( int in )
        {
            fnmscw = in;
        }
        void setMSCW( double iMSWC )
        {
            fmscw = iMSWC;
        }
        void setMSCWT( int i, double iMSWC, float iMSWC_T = -99. )
        {
            ftmscw[i] = iMSWC;
            ftmscw_sigma[i] = iMSWC_T;
        }
        void setMSCT( double iMSWT )
        {
            fmsct = iMSWT;
        }
        void setMSCTT( int i, double iMSTC, float iMSTC_T = -99. )
        {
            ftmsct[i] = iMSTC;
            ftmsct_sigma[i] = iMSTC_T;
        }
        void setMWR( double iM )
        {
            fmwr = iM;
        }
        void setMLR( double iM )
        {
            fmlr = iM;
        }
        void setMCDistanceToCameraCenter( double iMin, double iMax )
        {
            fMC_distance_to_cameracenter_min = iMin;
            fMC_distance_to_cameracenter_max = iMax;
        }
        void setDebug( bool iD )
        {
            fDebug = iD;
        }
        bool setOutputFile( string, string, string );
        bool useDispEnergy()
        {
            if( fDispAnalyzerEnergy )
            {
                return true;
            }
            return false;
        }
        bool terminate( TNamed* );                //!< write everything to disk
};
#endif
