/*! \class VTableLookupDataHandler
  \brief data class for mscw and energy reconstruction

  reads evndisp output trees, fill results from mscw and energy reconstruction

*/

#include "VTableLookupDataHandler.h"

VTableLookupDataHandler::VTableLookupDataHandler( bool iwrite, VTableLookupRunParameter* iT )
{
    fTLRunParameter = iT;
    if( !fTLRunParameter )
    {
        cout << "VTableLookupDataHandler::VTableLookupDataHandler error: to table lookup run parameters" << endl;
        exit( EXIT_FAILURE );
    }
    fDebug = fTLRunParameter->fDebug;
    if( fDebug > 0 )
    {
        cout << "VTableLookupDataHandler::VTableLookupDataHandler" << endl;
    }

    fwrite = iwrite;
    fNEntries = fTLRunParameter->fNentries;
    fEventDisplayFileFormat = 2;
    fTshowerpars = 0;
    fTshowerpars_QCCut = 0;
    fshowerpars = 0;
    fDeepLearnerpars = 0;
    fOTree = 0;
    fShortTree = fTLRunParameter->bShortTree;
    bWriteMCPars = fTLRunParameter->bWriteMCPars;
    fTreeWithParameterErrors = false;
    fNTel = 0;
    fNTelComb = 0;
    fTtelconfig = 0;
    foutputfile = "";

    fEmissionHeightCalculator = new VEmissionHeightCalculator();

    fEventStatus = true;

    // random number generator is needed only for random selection of events (optional)
    fRandom = new TRandom3();
    setSelectRandom( fTLRunParameter->fSelectRandom, fTLRunParameter->fSelectRandomSeed );

    fOutFile = 0;

    fNMethods = 0;
    fMethod = fTLRunParameter->rec_method;

    fIsMC = false;
    fMCEnergy = 0.;
    fZe = 0.;
    fEventCounter = 0;

    fEventWeight = 1.;
    // (hardwired DL parameters)
    // --> set true to read DL parameters
    // and fill them into the output tree
    fIsDeepLearner = false;

    /////////////////////////////////////////////////////////////////////////////////////
    // values needed by the optional stereo reconstruction
    fSSR_AxesAngles_min = 0.;
    fSSR_NImages_min    = 0;

    /////////////////////////////////////////////////////////////////////////////////////
    // weighting of energy spectrum
    // MC input spectrum
    fMCSpectralIndex = 2.0;
    fMCMinEnergy = 0.05;
    fMCMaxEnergy = 50.;
    // spectral index events are weighted to
    fSpectralIndex = 2.0;
    /////////////////////////////////////////////////////////////////////////////////////

    // MC spectra histograms
    hisList = new TList();
    hE0mc = new TH1D( "hE0mc", "MC energy spectrum", 1175, -2., 2.7 );
    hE0mc->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    hE0mc->SetYTitle( "number of events" );
    hE0mc->Sumw2();
    hisList->Add( hE0mc );

    hDE0mc = new TH2D( "hDE0mc", "distance vs. MC primary energy", 1175, -2., 2.7, 1000, 0., 2000. );
    hDE0mc->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    hDE0mc->SetYTitle( "distance to shower core [m]" );
    hDE0mc->SetZTitle( "number of events" );
    hDE0mc->Sumw2();
    hisList->Add( hDE0mc );

    hXYmc = new TH2D( "hXYmc", "MC core distribution", 2000, -2000., 2000., 2000, -2000., 2000. );
    hXYmc->Sumw2();
    hXYmc->SetXTitle( "x [m]" );
    hXYmc->SetYTitle( "y [m]" );
    hisList->Add( hXYmc );

    hWE0mc = new TH2D( "hWE0mc", "ang. dist. vs. energy", 1175, -2., 2.75, 500, 0., 10. );
    hWE0mc->SetTitle( "ang. distance vs. primary energy" );
    hWE0mc->SetXTitle( "log_{10} E_{MC} [TeV]" );
    hWE0mc->SetYTitle( "distance to camera center [deg]" );
    hWE0mc->SetZTitle( "number of showers" );
    hisList->Add( hWE0mc );

    hZe = new TH1D( "hZe", "cos(Ze)", 100, 0., 1.1 );
    hZe->SetXTitle( "cos(Ze)" );
    hisList->Add( hZe );

    // same with triggered events
    hE0trig = new TH1D( "hE0trig", "MC energy spectrum (triggered events)", 1175, -2., 2.7 );
    hE0trig->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    hE0trig->SetYTitle( "number of events" );
    hE0trig->Sumw2();
    hisList->Add( hE0trig );

    hDE0trig = new TH2D( "hDE0trig", "distance vs. MC primary energy (triggered events)", 1175, -2., 2.7, 1000, 0., 2000. );
    hDE0trig->SetXTitle( "log_{10} energy_{MC} [TeV]" );
    hDE0trig->SetYTitle( "distance to shower core [m]" );
    hDE0trig->SetZTitle( "number of events" );
    hDE0trig->Sumw2();
    hisList->Add( hDE0trig );

    hXYtrig = new TH2D( "hXYtrig", "core distribution (triggered events)", 2000, -2000., 2000., 2000, -2000., 2000. );
    hXYtrig->Sumw2();
    hXYtrig->SetXTitle( "x [m]" );
    hXYtrig->SetYTitle( "y [m]" );
    hisList->Add( hXYtrig );

    hWE0trig = new TH2D( "hWE0trig", "ang. dist. vs. energy", 1000, -2., 2, 500, 0., 10. );
    hWE0trig->SetTitle( "ang. distance vs. primary energy (triggered events)" );
    hWE0trig->SetXTitle( "log_{10} E_{MC} [TeV]" );
    hWE0trig->SetYTitle( "distance to camera center [deg]" );
    hWE0trig->SetZTitle( "number of showers" );
    hisList->Add( hWE0trig );

    // time cuts
    fMaxTotalTime = fTLRunParameter->fMaxRunTime;
    fTotalTime = 0.;
    fTotalTime0 = 0.;

    resetAll();

    fDispAnalyzerDirection = 0;
    fDispAnalyzerDirectionError = 0;
    fDispAnalyzerEnergy    = 0;
    fDispAnalyzerCore      = 0;
}

/*
 * fill results of analysis into output tree
 * (called data in the mscw file)
 */
void VTableLookupDataHandler::fill()
{
    if( !fOTree )
    {
        return;
    }

    // fill short variables
    unsigned int ii = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fImgSel_list[i] )
        {
            fR_short[ii]  = fR[i];
            fR_short_MC[ii]  = fR_MC[i];
            fES_short[ii] = fES[i];

            // cross for dispBDT debugging
            fcross_short[ii] = sqrt( ( fcen_y[i] + fYoff_intersect ) * ( fcen_y[i] + fYoff_intersect )
                                     + ( fcen_x[i] - fXoff_intersect ) * ( fcen_x[i] - fXoff_intersect ) );
            fcrossO_short[ii] = sqrt( ( fcen_y[i] + fYoff_edisp ) * ( fcen_y[i] + fYoff_edisp )
                                      + ( fcen_x[i] - fXoff_edisp ) * ( fcen_x[i] - fXoff_edisp ) );
            fdist_short[ii]   = fdist[i];
            fwidth_short[ii]  = fwidth[i];
            fdwidth_short[ii]  = fdwidth[i];
            flength_short[ii] = flength[i];
            fdlength_short[ii] = fdlength[i];
            fcen_x_short[ii] = fcen_x[i];
            fdcen_x_short[ii] = fdcen_x[i];
            fcen_y_short[ii] = fcen_y[i];
            fdcen_y_short[ii] = fdcen_y[i];
            fcosphi_short[ii] = fcosphi[i];
            fsinphi_short[ii] = fsinphi[i];
            fdphi_short[ii] = fdphi[i];
            floss_short[ii]   = floss[i];
            ffui_short[ii]    = ffui[i];
            fsize_short[ii]   = fsize[i];
            fntubes_short[ii] = fntubes[i];
            ftgrad_x_short[ii] = ftgrad_x[i];
            fasym_short[ii]  = fasym[i];
            fFitstat_short[ii] = fFitstat[i];

            ii++;
        }
    }

    if( fTLRunParameter->bWriteReconstructedEventsOnly >= 0
            || fTLRunParameter->bWriteReconstructedEventsOnly == -2 )
    {
        if( isReconstructed() )
        {
            fOTree->Fill();
        }
    }
    else
    {
        fOTree->Fill();
    }
}


void VTableLookupDataHandler::fillMChistograms()
{
    if( fIsMC )
    {
        // fill histograms with all simulated events
        double ilogE = log10( fMCEnergy );
        double idist = sqrt( fMCxcore * fMCxcore + fMCycore * fMCycore );
        double ioff  = sqrt( fMCxoff * fMCxoff + fMCyoff * fMCyoff );

        hE0mc->Fill( ilogE );
        hDE0mc->Fill( ilogE, idist );
        hWE0mc->Fill( ilogE, ioff );
        hXYmc->Fill( fMCxcore, fMCycore );
        hZe->Fill( fMCze );
        // fill histograms with all triggered events (require array trigger, at least 2 telescopes)
        if( fNTrig >= 2 )
        {
            hE0trig->Fill( ilogE );
            hDE0trig->Fill( ilogE, idist );
            hWE0trig->Fill( ilogE, ioff );
            hXYtrig->Fill( fMCxcore, fMCycore );
        }
    }
}


double VTableLookupDataHandler::getMCDistance()
{
    return sqrt( fMCxcore * fMCxcore + fMCycore * fMCycore );
}


/*!
    return values:

    true:     get successfully an event
    false:    end of data chain or time limit exceeded
*/
bool VTableLookupDataHandler::getNextEvent( bool bShort )
{
    if( fEventCounter < fNEntries && fTotalTime < fMaxTotalTime )
    {
        if( !randomSelected() )
        {
            fEventCounter++;
            return true;
        }
        fEventWeight = 1.;

        int next_event_status = 1;
        if( fEventDisplayFileFormat >= 2 )
        {
            next_event_status = fillNextEvent( bShort );
        }
        else
        {
            cout << "unknown eventdisplay file format: " << fEventDisplayFileFormat << endl;
            cout << "(possible old format? Format version: " << fEventDisplayFileFormat << ")" << endl;
            cout << "...exiting" << endl;
            exit( EXIT_FAILURE );
        }
        if( next_event_status == -1 )
        {
            return false;
        }

        // return false for non-valid (maybe not reconstructed?) event
        if( next_event_status == 0 )
        {
            return true;
        }

        // calculate theta2
        if( !fIsMC )
        {
            ftheta2 = ( fYoff_derot - fWobbleN ) * ( fYoff_derot - fWobbleN )
                      + ( fXoff_derot - fWobbleE ) * ( fXoff_derot - fWobbleE );
        }
        else
        {
            ftheta2 = ( fXoff - fMCxoff ) * ( fXoff - fMCxoff )
                      + ( fYoff - fMCyoff ) * ( fYoff - fMCyoff );
        }

        setEventWeightfromMCSpectrum();
    }
    else
    {
        return false;
    }
    return true;
}

/*
 * get next event from trees,
 * do quick reconstruction quality test,
 * fill variables
 * calculate missing variables
 *
 * returns -1 if no next event is found
 *
 */
int VTableLookupDataHandler::fillNextEvent( bool bShort )
{
    ///////////////////////////////////////////////////////////////////////////////
    // read partial event for quick reconstruction quality assessment
    if( !fTshowerpars_QCCut->GetEntry( fEventCounter ) )
    {
        return -1;
    }

    // count all events
    fNStats_All++;
    ////////////////////////////////////////////////////
    // read first all entries needed for run modes (filling and reading)
    if( fMethod >= ( int )fNMethods )
    {
        cout << "VTableLookupDataHandler::fillNextEvent() error, invalid array reconstruction record" << endl;
        cout << "\t maximum number of records are " << fNMethods << " (request is " << fMethod << ")" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    fNImages = fNImages_QCTree[fMethod];
    fchi2    = fchi2_QCTree[fMethod];

    // for table filling: check as soon as possible if the event is useful
    // (also expect that this will not change if stereo reconstruction
    //  is repeated; this is probably not true)
    if( fwrite && !isReconstructed( true ) )
    {
        fEventStatus = false;
        fEventCounter++;
        return 0;
    }
    ///////////////////////////////////
    // now read full event
    if( !fshowerpars->GetEntry( fEventCounter ) )
    {
        return -1;
    }
    if( fDeepLearnerpars )
    {
        // showerpars and deep learner trees should be in
        // sync; if not; feel default values
        if( !fDeepLearnerpars->GetEntry( fEventCounter ) )
        {
            dl_gammaness = -999.;
            dl_isGamma = false;
        }
    }
    // fill MC parameters
    if( fIsMC )
    {
        fMCEnergy = fshowerpars->MCe0;
        fMCaz     = fshowerpars->MCaz;
        fMCze     = fshowerpars->MCze;
        fMCxcore  = fshowerpars->MCxcore;
        fMCycore  = fshowerpars->MCycore;
        fMCxoff   = fshowerpars->MCxoff;
        fMCyoff   = fshowerpars->MCyoff;
        if( !bShort && !fShortTree && !fwrite )
        {
            fMCxcore_SC = fshowerpars->MCxcore_SC;
            fMCycore_SC = fshowerpars->MCycore_SC;
            fMCxcos = fshowerpars->MCxcos;
            fMCycos = fshowerpars->MCycos;
        }
        if( !fwrite )
        {
            fMCCorsikaRunID = fshowerpars->MCCorsikaRunID;
            fMCCorsikaShowerID = fshowerpars->MCCorsikaShowerID;
            fMCFirstInteractionHeight = fshowerpars->MCFirstInteractionHeight;
            fMCFirstInteractionDepth = fshowerpars->MCFirstInteractionDepth;
        }
    }
    fArrayPointing_Elevation = fshowerpars->ArrayPointing_Elevation;
    fArrayPointing_Azimuth   = fshowerpars->ArrayPointing_Azimuth;
    fArrayPointing_RotationAngle = fshowerpars->ArrayPointing_deRotationAngle_deg * TMath::DegToRad();

    // the following variables are not set in table filling mode
    if( !fwrite )
    {
        runNumber   = fshowerpars->runNumber;
        eventNumber = fshowerpars->eventNumber;
        if( fDebug > 1 )
        {
            cout << "===============================================================================" << endl;
            cout << "SHOWERPARS EVENT " << fshowerpars->eventNumber << "\t" << fEventCounter << "\t";
            cout << fshowerpars->NImages[fMethod] << "\t" << fshowerpars->Chi2[fMethod] << endl;
        }
        time = fshowerpars->Time;
        if( fEventCounter == 0 )
        {
            fTotalTime0 = time;
        }
        fTotalTime = time - fTotalTime0;

        for( unsigned int i = 0; i < fNTel; i++ )
        {
            fTelElevation[i] = fshowerpars->TelElevation[i];
            fTelAzimuth[i]   = fshowerpars->TelAzimuth[i];
        }
        if( !fIsMC )
        {
            MJD = fshowerpars->MJD;
        }
        fNTrig = 0;
        // determine number of triggered telescopes
        for( unsigned t = 0; t < fshowerpars->NTrig; t++)
        {
            unsigned int tel_trig = fshowerpars->Trig_list[t];
            if( t < fTLRunParameter->fTelToAnalyzeData.size() && fTLRunParameter->fTelToAnalyzeData[tel_trig] && fTLRunParameter->fTelToAnalyzeData[tel_trig]->fTelToAnalyze )
            {
                fNTrig++;
            }
        }

        if( !bShort )
        {
            LTrig = ( ULong64_t )fshowerpars->LTrig;
        }
        fDispDiff = fshowerpars->DispDiff[fMethod];
        fimg2_ang = fshowerpars->img2_ang[fMethod];
        if( !bShort && !fShortTree )
        {
            fRA = fshowerpars->ra[fMethod];
            fDec = fshowerpars->dec[fMethod];
            fstdS = fshowerpars->stds[fMethod];
            fXcore_SC = fshowerpars->Xcore_SC[fMethod];
            fYcore_SC = fshowerpars->Ycore_SC[fMethod];
            fstdP = fshowerpars->stdp[fMethod];
        }
    } // end (!fwrite)

    fZe = fshowerpars->Ze[fMethod];
    fAz = fshowerpars->Az[fMethod];
    fXcore = fshowerpars->Xcore[fMethod];
    fYcore = fshowerpars->Ycore[fMethod];
    // return if stereo reconstruction was not successful
    // (don't do this if stereo reconstruction is
    //  repeated)
    if( !fTLRunParameter->fRerunStereoReconstruction
            && ( TMath::IsNaN( fXcore ) || TMath::IsNaN( fYcore ) ) )
    {
        fXcore =  -999999.;
        fYcore =  -999999.;
        fEventCounter++;
        fEventStatus = false;
        if( fDebug > 1 )
        {
            cout << "\t RECONSTRUCTED CORE NAN" << endl;
        }
        return 0;
    }
    // standard stereo reconstruction
    fXoff = fshowerpars->Xoff[fMethod];
    fYoff = fshowerpars->Yoff[fMethod];
    fXoff_derot = fshowerpars->XoffDeRot[fMethod];
    fYoff_derot = fshowerpars->YoffDeRot[fMethod];

    ///////////////////////////////////////////////////////
    // set telescope selection variables

    // bit coded image selection
    // (note limitation in number of telescopes (<64)
    fImgSel = ( ULong64_t )fshowerpars->ImgSel[fMethod];
    unsigned int ii = 0;
    for( unsigned int i = 0; i < getNTelTypes(); i++ )
    {
        NImages_Ttype[i] = 0;
    }
    // list of selected telescopes
    // (but loop over all telescopes!)
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        fImgSel_list[i] = ( bool )fshowerpars->ImgSel_list[fMethod][i];

        // telescope removed via external list
        // e.g. using -sub_array_sim_telarray_counting
        if( fTLRunParameter && fTLRunParameter->fTelToAnalyzeData.size() == getNTel() )
        {
            if( fTLRunParameter->fTelToAnalyzeData[i] )
            {
                fweight[i]      = fTLRunParameter->fTelToAnalyzeData[i]->fWeight;
                fImgSel_list[i] = ( fImgSel_list[i] && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze );
            }
            else
            {
                fweight[i]      = 1.;
            }
            fweightDispBDTs[i] = fweight[i];   // set later
        }

        if( fImgSel_list[i] )
        {
            fImgSel_list_short[ii] = i;
            // count the number of telescopes of this type
            NImages_Ttype[getTelType_arraycounter( i )]++;
            ii++;
        }
    }
    // new total number of images
    fNImages = ( int )ii;

    // for filling of lookup tables: first do quality cuts, if not return
    if( fwrite )
    {
        if( !cut( true ) )
        {
            fEventCounter++;
            fEventStatus = false;
            if( fDebug > 1 )
            {
                cout << "\t CUT FAILED" << endl;
            }
            return 0;
        }
        else
        {
            fEventStatus = true;
        }
    }
    // (end of accessing showerpars tree)
    //////////////////////////////////////////

    ////////////////////////////////////////////
    // initialize tpars trees
    // loop over all telescopes
    bitset<8 * sizeof( unsigned long )> i_nimage; // for imagepattern
    i_nimage.reset();

    Double_t SizeFirstMax_temp = -1000.;
    Double_t SizeSecondMax_temp = -100.;
    ii = 0;
    unsigned int i_pixel_id0 = 0;
    PixelListNPixelNN = 0;
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        bool fReadTPars = false;
        if( i < ftpars.size() && ftpars[i] )
        {
            fReadTPars = true;
        }
        // check if the tpars for this telescope should be
        // read
        if( ( fTLRunParameter->bWriteReconstructedEventsOnly >= 0 )
                || fTLRunParameter->bWriteReconstructedEventsOnly == -2 || fwrite )
        {
            if( fImgSel_list[i] )
            {
                fReadTPars = true;
            }
            else
            {
                fReadTPars = false;
            }
        }
        // read only those telescope which were part of the reconstruction
        if( fReadTPars )
        {
            if( !ftpars[i] )
            {
                cout << "VTableLookupDataHandler::fillNextEvent error:";
                cout << "tree tpars not found (telescope " << i + 1 << ")" << endl;
                cout << "\t(run " << runNumber << ", " << eventNumber << ")" << endl;
                exit( EXIT_FAILURE );
            }
            ftpars[i]->GetEntry( fEventCounter );

            fdist[i] = ftpars[i]->dist;
            ffui[i] = ftpars[i]->fui;
            fsize[i] = ftpars[i]->size;
            fsize2[i] = ftpars[i]->size2;
            floss[i] = ftpars[i]->loss;
            ffracLow[i] = ftpars[i]->fracLow;
            fwidth[i] = ftpars[i]->width;
            flength[i] = ftpars[i]->length;
            ftgrad_x[i] = ftpars[i]->tgrad_x;
            if( i < fpointingCorrections.size() && fpointingCorrections[i]
                    && fpointingCorrections[i]->is_initialized() )
            {
                fpointingCorrections[i]->getEntry( fEventCounter );
                fcen_x[i] = fpointingCorrections[i]->getCorrected_cen_x( ftpars[i]->cen_x );
                fcen_y[i] = fpointingCorrections[i]->getCorrected_cen_y( ftpars[i]->cen_y );
                float phi = fpointingCorrections[i]->getCorrected_phi(
                                ftpars[i]->cen_x,
                                ftpars[i]->cen_y,
                                ftpars[i]->f_d,
                                ftpars[i]->f_s,
                                ftpars[i]->f_sdevxy );
                fcosphi[i] = cos( phi );
                fsinphi[i] = sin( phi );
            }
            else
            {
                fcen_x[i] = ftpars[i]->cen_x;
                fcen_y[i] = ftpars[i]->cen_y;
                fcosphi[i] = ftpars[i]->cosphi;
                fsinphi[i] = ftpars[i]->sinphi;
            }
            if( ftpars[i]->hasParameterErrors() )
            {
                fTreeWithParameterErrors = true;
                fdcen_x[i] = ftpars[i]->dcen_x;
                fdcen_y[i] = ftpars[i]->dcen_y;
                fdlength[i] = ftpars[i]->dlength;
                fdwidth[i] = ftpars[i]->dwidth;
                fdphi[i] = ftpars[i]->dphi;
            }
            if( ftpars[i]->hasPixelList() && fTLRunParameter->fWritePixelLists )
            {
                PixelListN[ii] = ftpars[i]->PixelListN;
                PixelListNPixelNN += ftpars[i]->PixelListN;
                unsigned int i_pix_id = 0;
                for( unsigned int p = 0; p < PixelListN[ii]; p++ )
                {
                    i_pix_id = i_pixel_id0 + p;
                    PixelID[i_pix_id] = ftpars[i]->PixelID[p];
                    PixelType[i_pix_id] = ftpars[i]->PixelType[p];
                    PixelIntensity[i_pix_id] = ftpars[i]->PixelIntensity[p];
                    PixelTimingT0[i_pix_id] = ftpars[i]->PixelTimingT0[p];
                    PixelPE[i_pix_id] = ftpars[i]->PixelPE[p];
                }
                i_pixel_id0 += PixelListN[ii];
                ii++;
            }

            if( fsize[i] > SizeSecondMax_temp )
            {
                if( fsize[i] > SizeFirstMax_temp )
                {
                    SizeSecondMax_temp = SizeFirstMax_temp;
                    SizeFirstMax_temp = fsize[i];
                }
                else
                {
                    SizeSecondMax_temp = fsize[i];
                }
            }
            fCurrentNoiseLevel[i] = ftpars[i]->meanPedvar_Image;
            if( !bShort )
            {
                fmeanPedvar_ImageT[i] = ftpars[i]->meanPedvar_Image;
                fntubes[i] = ftpars[i]->ntubes;
                fnsat[i] = ftpars[i]->nsat;
                fnlowgain[i] = ftpars[i]->nlowgain;
                falpha[i] = ftpars[i]->alpha;
                flos[i] = ftpars[i]->los;
                fasym[i] = ftpars[i]->asymmetry;
                fmax1[i] = ftpars[i]->max[0];
                fmax2[i] = ftpars[i]->max[1];
                fmax3[i] = ftpars[i]->max[2];
                fmaxindex1[i] = ftpars[i]->index_of_max[0];
                fmaxindex2[i] = ftpars[i]->index_of_max[1];
                fmaxindex3[i] = ftpars[i]->index_of_max[2];
                ftchisq_x[i] = ftpars[i]->tchisq_x;
                fFitstat[i] = ftpars[i]->Fitstat;
            }
        }
        else
        {
            if( !fwrite )
            {
                resetImageParameters( i );
            }
        }

        // bit coding for telescope used in analysis
        // (small arrays only)
        if( !bShort && fntubes[i] > 4 && i < i_nimage.size() && i < 10 )
        {
            i_nimage.set( i, 1 );
        }
    }
    fmeanPedvar_Image = calculateMeanNoiseLevel( true );

    if( SizeSecondMax_temp > 0. )
    {
        fSizeSecondMax = SizeSecondMax_temp;
    }
    if( fNImages == 1 )
    {
        fSizeSecondMax =  SizeFirstMax_temp;
    }

    ///////////////////////////////////////////////////////////
    // calculate distances
    calcDistances();

    ///////////////////////////////////////////////////////////
    // calculate emission height (not for writing of tables)
    if( fNImages > 1 && !fwrite )
    {
        calcEmissionHeights();
    }
    else
    {
        fEmissionHeightMean = 1.e-10;
        fEmissionHeightChi2 = 1.e-10;
        fNTelPairs = 0;
    }


    //////////////////////////////////////////////////////////
    // !!! SPECIAL AND EXPERT USAGE ONLY !!!
    // redo the stereo (direction and core) reconstruction
    // Works for MC only!!
    //
    if( fTLRunParameter->fRerunStereoReconstruction )
    {
        doStereoReconstruction();
    }

    //////////////////////////////////////////////////////////
    // !!! SPECIAL AND EXPERT USAGE ONLY !!!
    // dispCore
    // core reconstruction using the disp MVA
    // This is preliminary and works for MC events only!
    //
    if( fDispAnalyzerCore )
    {
        fDispAnalyzerCore->setQualityCuts( fSSR_NImages_min, fSSR_AxesAngles_min,
                                           fTLRunParameter->fmaxdist,
                                           fTLRunParameter->fmaxloss,
                                           fTLRunParameter->fminfui,
                                           fmaxdist_qc );
        fDispAnalyzerCore->calculateCore(
            getNTel(),
            fArrayPointing_Elevation, fArrayPointing_Azimuth,
            fTelX, fTelY, fTelZ,
            fTel_type,
            getSize( 1., true, false ),
            fcen_x, fcen_y,
            fcosphi, fsinphi,
            fwidth, flength,
            fasym, ftgrad_x,
            floss, fntubes,
            getWeight(),
            fXoff, fYoff,
            getDistanceToCore(),
            fXcore, fYcore,
            fXoff, fYoff,
            ffui );

        // fill results
        for( unsigned int i = 0; i < getNTel(); i++ )
        {
            fR[i] = fDispAnalyzerCore->getCoreDistance( i );
            // (strictly not correct, but good enough)
            fRTel[i] = fDispAnalyzerCore->getCoreDistance( i );
        }
    }

    //////////////////////////////////////////////////////////
    // !!! SPECIAL AND EXPERT USAGE ONLY !!!
    // dispEnergy
    // energy reconstruction using the disp MVA
    // This is preliminary and works for MC events only!
    //
    if( fDispAnalyzerEnergy )
    {
        fDispAnalyzerEnergy->setQualityCuts( fSSR_NImages_min, fSSR_AxesAngles_min,
                                             fTLRunParameter->fmaxdist,
                                             fTLRunParameter->fmaxloss,
                                             fTLRunParameter->fminfui,
                                             fmaxdist_qc );
        fDispAnalyzerEnergy->calculateEnergies(
            getNTel(),
            fArrayPointing_Elevation, fArrayPointing_Azimuth,
            fTel_type,
            getSize( 1., true, false ),
            fcen_x, fcen_y,
            fcosphi, fsinphi,
            fwidth, flength,
            fasym, ftgrad_x,
            floss, fntubes,
            getWeight(true),
            fXoff, fYoff,
            getDistanceToCoreTel(),
            fEmissionHeightMean,
            fMCEnergy,
            ffui );

        // fill results
        setEnergy( fDispAnalyzerEnergy->getEnergy(),
                   fDispAnalyzerEnergy->getEnergyChi2(),
                   fDispAnalyzerEnergy->getEnergydES(),
                   fDispAnalyzerEnergy->getEnergyMedianAbsoluteError() );

        for( unsigned int i = 0; i < getNTel(); i++ )
        {
            setEnergyT( i, fDispAnalyzerEnergy->getEnergyT( i ) );
        }
        setNEnergyT( fDispAnalyzerEnergy->getEnergyNT() );
        setNEnergyQuality( fDispAnalyzerEnergy->getEnergyQualityLabel() );
    }

    fEventCounter++;
    return 1;
}

/*
 * redo stereo reconstruction (core and direction)
 *
 * this works for MC only
 * not all stereo reconstruction methods are implemented
 * (quick and dirty implementation for CTA)
 *
 * does not take into account pointing corrections
 * (as e.g. given by the VPM)
*/
void VTableLookupDataHandler::doStereoReconstruction()
{
    // save original values
    fXoff_edisp = fXoff;
    fYoff_edisp = fYoff;
    ///////////////////////////
    // stereo reconstruction
    // (rcs_method4)
    VSimpleStereoReconstructor i_SR;
    // minimal value; just used to initialize disp method
    i_SR.initialize( fSSR_NImages_min, fSSR_AxesAngles_min );
    i_SR.reconstruct_direction_and_core( getNTel(),
                                         fArrayPointing_Elevation, fArrayPointing_Azimuth,
                                         fTelX, fTelY, fTelZ,
                                         getSize( 1., true, false ),
                                         fcen_x, fcen_y,
                                         fcosphi, fsinphi,
                                         fwidth, flength,
                                         getWeight() );
    // store results from line intersection for debugging
    fXoff_intersect = i_SR.fShower_Xoffset;
    fYoff_intersect = i_SR.fShower_Yoffset;

    ////////////////////////////////////////////////////////////////////
    // DISP method for updated disp reconstruction
    ////////////////////////////////////////////////////////////////////
    if( fDispAnalyzerDirection
            && fNImages <= ( int )fTLRunParameter->fRerunStereoReconstruction_BDTNImages_max )
    {

        vector< float > iDispError( getNTel(), -9999. );

        ////////////////////////////////////////////////////////////////////
        // estimate error on direction reconstruction from DISP method
        ////////////////////////////////////////////////////////////////////
        if( fDispAnalyzerDirectionError )
        {
            fDispAnalyzerDirectionError->calculateExpectedDirectionError(
                getNTel(),
                fArrayPointing_Elevation, fArrayPointing_Azimuth,
                fTel_type,
                getSize( 1., true, false ),
                fcen_x, fcen_y,
                fcosphi, fsinphi,
                fwidth, flength,
                fasym, ftgrad_x,
                floss, fntubes,
                getWeight(),
                i_SR.fShower_Xoffset, i_SR.fShower_Yoffset,
                ffui );

            // get estimated error on direction reconstruction
            for( unsigned int t = 0; t < getNTel(); t++ )
            {
                iDispError[t] = fDispAnalyzerDirectionError->getDispErrorT( t );
                fweightDispBDTs[t] = iDispError[t];
            }
        }

        // use weighting calculated from disp error
        fDispAnalyzerDirection->setDispErrorWeighting( fDispAnalyzerDirectionError != 0,
                fTLRunParameter->fDispError_BDTWeight );
        fDispAnalyzerDirection->setQualityCuts( fSSR_NImages_min, fSSR_AxesAngles_min,
                                                fTLRunParameter->fmaxdist,
                                                fTLRunParameter->fmaxloss,
                                                fTLRunParameter->fminfui,
                                                fmaxdist_qc );
        fDispAnalyzerDirection->calculateMeanDirection(
            getNTel(),
            fArrayPointing_Elevation, fArrayPointing_Azimuth,
            fTel_type,
            getSize( 1., true, false ),
            fcen_x, fcen_y,
            fcosphi, fsinphi,
            fwidth, flength,
            fasym, ftgrad_x,
            floss, fntubes,
            getWeight(),
            i_SR.fShower_Xoffset, i_SR.fShower_Yoffset,
            iDispError, ffui );
        // reconstructed direction by disp method:
        fimg2_ang = fDispAnalyzerDirection->getAngDiff();
        fnxyoff = fDispAnalyzerDirection->getXYWeight_disp().size();
        if( fnxyoff > 0 )
        {
            fXoff = fDispAnalyzerDirection->getXcoordinate_disp();
            fYoff = fDispAnalyzerDirection->getYcoordinate_disp();

            // dispersion of disp values
            fDispDiff = fDispAnalyzerDirection->getDispDiff();
            fchi2 = fDispDiff;
            // for az / ze calculation
            i_SR.fillShowerDirection( fXoff, fYoff );
            for( unsigned int t = 0; t < fnxyoff; t++ )
            {
                fXoff_T[t] = fDispAnalyzerDirection->getXcoordinate_disp( t );
                fYoff_T[t] = fDispAnalyzerDirection->getYcoordinate_disp( t );
                fWoff_T[t] = fDispAnalyzerDirection->getXYWeight_disp( t );
                fDoff_T[t] = fDispAnalyzerDirection->get_disp( t );
                fToff_T[t] = fDispAnalyzerDirection->get_disp_tel_list( t );
            }
        }
        else
        {
            fXoff = -99.;
            fYoff = -99.;
        }
        fDispAbsSumWeigth = fDispAnalyzerDirection->get_abs_sum_disp_weight();
    }
    ////////////////////////////////////////////////////////////////////
    // Standard (intersection) method for all other cases
    ////////////////////////////////////////////////////////////////////
    else
    {
        fXoff  = i_SR.fShower_Xoffset;
        fYoff  = i_SR.fShower_Yoffset;
        fstdS  = i_SR.fShower_stdS;
        fchi2  = i_SR.fShower_Chi2;
        fDispDiff = i_SR.fShower_DispDiff;
        fimg2_ang = i_SR.fiangdiff;
        fXoff_derot = i_SR.fShower_Xoffset;
        fYoff_derot = i_SR.fShower_Yoffset;
        fstdS = i_SR.fShower_DispDiff;
        fDispAbsSumWeigth = 0.;
    }

    // overwrite the values read from the evndisp file with the newly
    // calculated values
    if( fIsMC )
    {
        fXoff_derot = fXoff; // MC only!
        fYoff_derot = fYoff; // MC only!
    }
    // derotate coordinates
    else
    {
        fXoff_derot = fXoff * cos( fArrayPointing_RotationAngle )
                      - fYoff * sin( fArrayPointing_RotationAngle );
        fYoff_derot = fYoff * cos( fArrayPointing_RotationAngle )
                      + fXoff * sin( fArrayPointing_RotationAngle );
    }
    fZe    = i_SR.fShower_Ze;
    fAz    = i_SR.fShower_Az;
    fXcore = i_SR.fShower_Xcore;
    fYcore = i_SR.fShower_Ycore;
}

/*
 * check input data / chains for consistency
 *
 * chains marked as 'recovered' by root cannot be used,
 * as usually the analysis does not complete correctly
 * for these chains
 *
*/
bool VTableLookupDataHandler::checkIfFilesInChainAreRecovered( TChain* c )
{
    if( !c )
    {
        cout << "VTableLookupDataHandler::checkIfFilesInChainAreRecovered() error: no chain" << endl;
        return true;
    }

    TObjArray* fileElements = c->GetListOfFiles();
    if( !fileElements )
    {
        cout << "VTableLookupDataHandler::checkIfFilesInChainAreRecovered() error: no files in chain" << endl;
        return true;
    }
    TChainElement* chEl = 0;
    TIter next( fileElements );
    while( ( chEl = ( TChainElement* )next() ) )
    {
        TFile* ifInput = new TFile( chEl->GetTitle() );
        if( ifInput->IsZombie() )
        {
            cout << "VTableLookupDataHandler::checkIfFilesInChainAreRecovered() error: file cannot be recovered; possibly not complete" << endl;
            cout << "\t " << chEl->GetTitle() << endl;
            return true;
        }
        if( ifInput->TestBit( TFile::kRecovered ) )
        {
            cout << "VTableLookupDataHandler::checkIfFilesInChainAreRecovered() error: file recovered; possibly not complete" << endl;
            cout << "\t " << chEl->GetTitle() << endl;
            return true;
        }
        ifInput->Close();
    }

    return false;
}

/*
 * print list of telescopes
 *
 * iUsedTelescopes = 0: print all telescopes
 * iUsedTelescopes = 1: print telescope used only
 * iUsedTelescopes = 2: print telescope not used only
 *
 */
void VTableLookupDataHandler::printTelescopesList( unsigned int iPrintParameter )
{
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        bool bUsed = true;
        if( fTLRunParameter && fTLRunParameter->fTelToAnalyzeData.size() == getNTel()
                && fTLRunParameter->fTelToAnalyzeData[i]
                && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
        {
            bUsed = true;
        }
        else
        {
            bUsed = false;
        }
        if( iPrintParameter == 1 && !bUsed )
        {
            continue;
        }
        if( iPrintParameter == 2 && bUsed )
        {
            continue;
        }
        // print list

        if( bUsed )
        {
            cout << "\t TELESCOPE " << i + 1 << "\t";
        }
        else
        {
            cout << "\t telescope " << i + 1 << "\t";
        }
        ios::fmtflags f( cout.flags() );
        std::streamsize ss = std::cout.precision();
        cout << setprecision( 2 ) << fixed;
        cout << "x:" << fTelX[i] << " [m]\ty:" << fTelY[i] << " [m]\tz:" << fTelZ[i] << " [m]\t";
        cout << setprecision( ss ) << std::resetiosflags( std::ios::showbase );
        cout << "type " << fTel_type[i];
        cout.flags( f );

        if( bUsed )
        {
            cout << " (type counter " << getTelType_arraycounter( i ) << ")";
            cout << ", used in analysis";
            if( fTLRunParameter->fRerunStereoReconstruction )
            {
                cout << " (stereo weight: " << fTLRunParameter->fTelToAnalyzeData[i]->fWeight;
                cout << ", dist cut: " << fmaxdist_qc[i];
                cout << ")";
            }
        }
        else
        {
            cout << ", not used in analysis";
        }
        cout << endl;
    }
}

/*
* setup input file chains for data
*
* read all telescope configuration and
* run parameters from disk
*
* set input and output trees
*
*/
bool VTableLookupDataHandler::setInputFile( vector< string > iInput )
{
    finputfile = iInput;
    // need to find suffix .root: add it if it doesn't exist
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        if( finputfile[i].find( ".root" ) == string::npos )
        {
            cout << "TableLookupDataHandler::setInputFile: adding .root suffix to file name" << endl;
            finputfile[i] += ".root";
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // get telescope configuration
    // get it from the telescope configuration tree (if available), else assume two telescope setup
    fTtelconfig = new TChain( "telconfig" );
    int iNFil_sum = 0;
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        int iNFil = fTtelconfig->Add( finputfile[i].c_str() );
        if( iNFil == 0 )
        {
            cout << "error: no file(s) in chain" << endl;
            exit( EXIT_FAILURE );
        }
        iNFil_sum += iNFil;
    }
    cout << iNFil_sum << " file(s) in chain " << endl;
    // don't check each file for CTA sims -> this is very inefficient and it takes a long time
    if( !fTLRunParameter->fPE && !fTLRunParameter->fWriteTables )
    {
        if( checkIfFilesInChainAreRecovered( fTtelconfig ) )
        {
            cout << "VTableLookupDataHandler::setInputFile() error: some file are not properly closed" << endl;
            cout << "exit..." << endl;
            exit( EXIT_FAILURE );
        }
    }

    ////////////////////////////////////
    // read in telescope configuration
    fList_of_Tel_type.clear();
    if( fTtelconfig )
    {
        ftelconfig = new Ctelconfig( fTtelconfig );
        ftelconfig->GetEntry( 0 );
        fNTel = ftelconfig->NTel;
        if( fNTel > getMaxNbrTel() )
        {
            cout << "VTableLookupDataHandler::setInputFile: error too many telescopes " << fNTel << "\t" << getMaxNbrTel() << endl;
            exit( EXIT_FAILURE );
        }
        fNTelComb = ( unsigned int )TMath::Nint( TMath::Power( 2., ( double )fNTel ) );
        for( unsigned int i = 0; i < fNTel; i++ )
        {
            ftelconfig->GetEntry( i );
            fTelX[i] = ftelconfig->TelX;
            fTelY[i] = ftelconfig->TelY;
            fTelZ[i] = ftelconfig->TelZ;
            fFocalLength[i] = ftelconfig->FocalLength;
            fTelFOV[i]      = ftelconfig->FOV;
            fTel_type[i]    = ftelconfig->TelType;
            // only count telescopes which are in the list of telescopes
            // used in the analysis
            if( fTLRunParameter && i < fTLRunParameter->fTelToAnalyzeData.size()
                    && fTLRunParameter->fTelToAnalyzeData[i] && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                if( fList_of_Tel_type.find( ftelconfig->TelType ) != fList_of_Tel_type.end() )
                {
                    fList_of_Tel_type[ftelconfig->TelType]++;
                }
                else
                {
                    fList_of_Tel_type[ftelconfig->TelType] = 1;
                }
            }
            else if( fTLRunParameter && fTLRunParameter->fTelToAnalyzeData.size() == 0 )
            {
                if( fList_of_Tel_type.find( ftelconfig->TelType ) != fList_of_Tel_type.end() )
                {
                    fList_of_Tel_type[ftelconfig->TelType]++;
                }
                else
                {
                    fList_of_Tel_type[ftelconfig->TelType] = 1;
                }
            }
            // quality cut on maximal distance to camera centre
            if( fTLRunParameter->fmaxdistfraction > 0. )
            {
                fmaxdist_qc[i] = fTLRunParameter->fmaxdistfraction * 0.5 * fTelFOV[i];
            }
            else
            {
                fmaxdist_qc[i] = fTLRunParameter->fmaxdist;
            }
        }
        // special case: no telescope given in teltoanalyze vector
        // number of different telescope types
        fNTelTypes = ( int )fList_of_Tel_type.size();
        // array with telescope type IDs
        unsigned int z = 0;
        for( fList_of_Tel_type_iterator = fList_of_Tel_type.begin();
                fList_of_Tel_type_iterator != fList_of_Tel_type.end();
                fList_of_Tel_type_iterator++ )
        {
            fTtypeID[z] = fList_of_Tel_type_iterator->first;
            z++;
        }
    }
    else
    {
        cout << "VTableLookupDataHandler::setInputFile error: no telescope configurations found " << endl;
        cout << "...exiting" << endl;
        exit( EXIT_FAILURE );
    }

    ////////////////////////////////////////////
    // print telescope configuration to the screen
    int i_ntel_on = fNTel;

    // count number of telescopes to be included in the analysis
    if( fTLRunParameter && fTLRunParameter->fTelToAnalyzeData.size() == getNTel() )
    {
        i_ntel_on = 0;
        for( unsigned int i = 0; i < fTLRunParameter->fTelToAnalyzeData.size(); i++ )
        {
            if( fTLRunParameter->fTelToAnalyzeData[i] && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                i_ntel_on++;
            }
        }
    }
    cout << endl << "total number of telescopes: " << fNTel;
    cout << " (" << i_ntel_on << " included in analysis)" << endl;
    initializeTelTypeVector();

    // print list of telescopes used
    cout << "Telescope used in analysis: " << endl;
    printTelescopesList( 1 );
    // print list of telescopes not used
    cout << endl << "Telescope not used in analysis: " << endl;
    printTelescopesList( 2 );
    cout << endl;
    cout << "list of telescope types (" << fList_of_Tel_type.size() << "): ";
    for( fList_of_Tel_type_iterator = fList_of_Tel_type.begin(); fList_of_Tel_type_iterator != fList_of_Tel_type.end();
            fList_of_Tel_type_iterator++ )
    {
        cout << "  " << fList_of_Tel_type_iterator->first << " (" << fList_of_Tel_type_iterator->second << " telescopes)";
    }
    cout << endl;
    // telconfig trees are not used anymore
    delete ftelconfig;
    ftelconfig = 0;
    // fTtelconfig = 0;
    // (end of telescope configuration)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // copy telescope positions to emission height calculator
    fEmissionHeightCalculator->setTelescopePositions( fNTel, fTelX, fTelY, fTelZ );

    // define trigger histogram
    char iName[100];
    char iDir[1000];
    unsigned int bShort = false;
    // get shower parameter tree
    fTshowerpars = new TChain( "showerpars" );
    fTshowerpars_QCCut = new TChain( "showerpars" );
    if( fIsDeepLearner )
    {
        fDeepLearnerpars = new TChain( "data_DL" );
    }

    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        fTshowerpars->Add( finputfile[i].c_str() );
        fTshowerpars_QCCut->Add( finputfile[i].c_str() );
        if( fDeepLearnerpars )
        {
            fDeepLearnerpars->Add( finputfile[i].c_str() );
        }
    }
    if( fDeepLearnerpars )
    {
        fDeepLearnerpars->SetBranchAddress( "dl_gammaness", &dl_gammaness );
        fDeepLearnerpars->SetBranchAddress( "dl_isGamma", &dl_isGamma );
    }
    if( !fTshowerpars )
    {
        cout << "VTableLookupDataHandler::setInputFile: error while retrieving data trees (2)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // check validity of showerpars tree
    if( !fTshowerpars->GetBranchStatus( "runNumber" ) )
    {
        cout << "VTableLookupDataHandler::setInputFile: error while retrieving data trees (2b)" << endl;
        cout << "exiting..." << endl;
        exit( EXIT_FAILURE );
    }
    // check if input data is MC
    if( fTshowerpars->GetBranchStatus( "MCe0" ) )
    {
        fIsMC = true;
        cout << "input data is of Monte Carlo type" << endl;
    }
    else
    {
        fIsMC = false;
        cout << "input data is not Monte Carlo type" << endl;
    }
    // for table filling: minimizing reading of trees
    if( fwrite )
    {
        bShort = 2;
    }

    // update runparameters
    fTLRunParameter->update( fTshowerpars );
    // get file format version of eventdisplay (tree version)
    if( fTLRunParameter )
    {
        fEventDisplayFileFormat = fTLRunParameter->getEVNDISP_TREE_VERSION();
        bShort                  = ( unsigned int )fTLRunParameter->getEVNDISP_TREE_isShort( fTshowerpars->GetTree() );
    }
    // check file format and initialize trees
    if( fEventDisplayFileFormat >= 2 )
    {
        if( bShort )
        {
            cout << "input data is of eventdisplay short tree output format (" << bShort << ")" << endl;
        }
        fshowerpars = new Cshowerpars( fTshowerpars, fIsMC, bShort );
        fIsMC = fshowerpars->isMC();
    }
    else
    {
        fEventDisplayFileFormat = 1;
    }
    // initialize quality cut chain
    fTshowerpars_QCCut->SetBranchAddress( "NMethods", &fNMethods );
    fTshowerpars_QCCut->SetBranchAddress( "NImages", fNImages_QCTree );
    fTshowerpars_QCCut->SetBranchAddress( "Chi2", fchi2_QCTree );

    // get individual image parameter trees
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        TChain* iT = new TChain( "tpars" );
        sprintf( iName, "pointing_%u", i + 1 );
        // pointing correction chain
        TChain* iPC = new TChain( iName );
        for( unsigned int f = 0; f < finputfile.size(); f++ )
        {
            sprintf( iDir, "%s/Tel_%u/tpars", finputfile[f].c_str(), i + 1 );
            iT->Add( iDir );
            // no pointing corrections for MC analysis
            if( !fIsMC )
            {
                sprintf( iDir, "%s/Tel_%u/pointing_%u", finputfile[f].c_str(), i + 1, i + 1 );
                iPC->Add( iDir );
            }
        }
        if( !iT )
        {
            cout << "VTableLookupDataHandler::setInputFile: error while retrieving data trees (3)" << endl;
            exit( EXIT_FAILURE );
        }
        // get first entry to check if chain is there
        // gErrorIgnoreLevel = 5000;
        if( iT->GetEntry( 0 ) > 0 )
        {
            if( fEventDisplayFileFormat >= 2 )
            {
                if( fEventDisplayFileFormat < 5 )
                {
                    if( iT->GetBranchStatus( "loss" ) )
                    {
                        fEventDisplayFileFormat = 3;
                    }
                    if( iT->GetBranchStatus( "meanPedvar_Image" ) )
                    {
                        fEventDisplayFileFormat = 5;
                    }
                }
                ftpars.push_back( new Ctpars( iT, fIsMC, bShort ) );
            }
        }
        else
        {
            ftpars.push_back( 0 );
        }
        fpointingCorrections.push_back( new VPointingCorrectionsTreeReader( iPC ) );
        gErrorIgnoreLevel = 0;
    }
    cout << "reading eventdisplay file format version " << fEventDisplayFileFormat;
    if( fIsMC )
    {
        cout << " (source files are Monte Carlo)";
    }
    cout << endl;

    ////////////////////////////////////////////////////////////////////////////////////
    // calculating median of pedvar distribution (not if input data is of PE format)
    fNoiseLevel.clear();
    fCurrentNoiseLevel.assign( fNTel, 0. );
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        // standard data format
        if( !fTLRunParameter->fPE )
        {
            if( fDebug > 1 )
            {
                cout << "VTableLookupDataHandler::setInputFile() calculating pedvar for telescope " << i + 1 << endl;
            }
            sprintf( iName, "calib_%u", i + 1 );
            TChain iPedVars( iName );
            for( unsigned int f = 0; f < finputfile.size(); f++ )
            {
                gErrorIgnoreLevel = 5000;
                sprintf( iDir, "%s/Tel_%u/calib_%u", finputfile[f].c_str(), i + 1, i + 1 );
                if( !iPedVars.Add( iDir ) )
                {
                    cout << "VTableLookupDataHandler::setInputFile: error while retrieving pedvars trees" << endl;
                    cout << "exiting..." << endl;
                    exit( EXIT_FAILURE );
                }
                if( iPedVars.GetEntries() == 0 )
                {
                    // backwards compatibility: read calibration tree from a different directory (note: this produces a root error message)
                    sprintf( iDir, "%s/Tel_%u/calibration/calib_%u", finputfile[f].c_str(), i + 1, i + 1 );
                    if( !iPedVars.Add( iDir ) )
                    {
                        cout << "VTableLookupDataHandler::setInputFile: error while retrieving pedvars trees" << endl;
                        cout << "exiting..." << endl;
                        exit( EXIT_FAILURE );
                    }
                }
                gErrorIgnoreLevel = 0;
            }
            gErrorIgnoreLevel = 5000;
            double pedvar = 0.;
            double mpedvar = 0.;
            double npedvar = 0.;
            double state = 0;
            iPedVars.SetBranchAddress( "pedvar", &pedvar );
            if( iPedVars.GetBranchStatus( "state" ) )
            {
                iPedVars.SetBranchAddress( "state", &state );
            }

            sprintf( iName, "ht_%u", i + 1 );
            TH1D h( iName, "", 1000, 0., 50. );

            if( fDebug > 1 )
            {
                cout << "VTableLookupDataHandler::setInputFile() calculating pedvar for telescope ";
                cout << i + 1 << ", number of entries: " << iPedVars.GetEntries() << endl;
            }
            for( int n = 0; n < iPedVars.GetEntries(); n++ )
            {
                iPedVars.GetEntry( n );

                if( pedvar > 0. && state == 0 )
                {
                    mpedvar += pedvar;
                    npedvar++;
                    h.Fill( pedvar );
                }
            }
            double xq[1];
            double yq[1];
            xq[0] = 0.5;
            yq[0] = 0.;
            if( h.GetEntries() > 0. )
            {
                h.GetQuantiles( 1, yq, xq );
            }
            fNoiseLevel.push_back( yq[0] );
            gErrorIgnoreLevel = 0;
            if( fDebug > 1 )
            {
                cout << "VTableLookupDataHandler::setInputFile() calculating pedvar for telescope (results): " << i + 1 << "\t" << yq[0] << endl;
            }
        }
        // PE format -> ignore noise level calculation
        else
        {
            fNoiseLevel.push_back( 0. );
        }
    }


    // temporary list of telescopes required for disp analysers
    vector<ULong64_t> i_TelTypeList;
    for( fList_of_Tel_type_iterator = fList_of_Tel_type.begin();
            fList_of_Tel_type_iterator != fList_of_Tel_type.end();
            fList_of_Tel_type_iterator++ )
    {
        if( fList_of_Tel_type_iterator->second > 0 )
        {
            i_TelTypeList.push_back( fList_of_Tel_type_iterator->first );
        }
    }

    /////////////////////////////////////////
    // initialize Disp Analyzer for direction reconstruction
    // (if required)
    if( fTLRunParameter->fRerunStereoReconstruction_BDTFileName.size() > 0. )
    {
        cout << endl;
        cout << "Initializing BDT disp analyzer for direction reconstruction" << endl;
        cout << "===========================================================" << endl << endl;
        fDispAnalyzerDirection = new VDispAnalyzer();
        fDispAnalyzerDirection->setTelescopeTypeList( i_TelTypeList );
        fDispAnalyzerDirection->initialize( fTLRunParameter->fRerunStereoReconstruction_BDTFileName, "TMVABDT" );
    }
    /////////////////////////////////////////
    // initialize Disp Analyzer for error on direction reconstruction
    // (if required)
    if( fTLRunParameter->fDispError_BDTFileName.size() > 0. )
    {
        cout << endl;
        cout << "Initializing BDT disp analyzer for estimation of disp error" << endl;
        cout << "===========================================================" << endl << endl;
        cout << "\t error weighting parameter: " << fTLRunParameter->fDispError_BDTWeight << endl;
        fDispAnalyzerDirectionError = new VDispAnalyzer();
        fDispAnalyzerDirectionError->setTelescopeTypeList( i_TelTypeList );
        fDispAnalyzerDirectionError->initialize( fTLRunParameter->fDispError_BDTFileName, "TMVABDT", "BDTDispError" );
    }
    /////////////////////////////////////////
    // initialize Disp Analyzer for core reconstruction
    // (if required)
    if( fTLRunParameter->fCoreReconstruction_BDTFileName.size() > 0. )
    {
        cout << endl;
        cout << "Initializing BDT disp analyzer for core reconstruction" << endl;
        cout << "===========================================================" << endl << endl;
        fDispAnalyzerCore = new VDispAnalyzer();
        fDispAnalyzerCore->setTelescopeTypeList( i_TelTypeList );
        fDispAnalyzerCore->initialize( fTLRunParameter->fCoreReconstruction_BDTFileName, "TMVABDT", "BDTDispCore" );
    }
    /////////////////////////////////////////
    // initialize Disp Analyzer for energy reconstruction
    // (if required)
    if( fTLRunParameter->fEnergyReconstruction_BDTFileName.size() > 0. )
    {
        cout << endl;
        cout << "Initializing BDT disp analyzer for energy reconstruction" << endl;
        cout << "===========================================================" << endl << endl;
        fDispAnalyzerEnergy = new VDispAnalyzer();
        fDispAnalyzerEnergy->setTelescopeTypeList( i_TelTypeList );
        fDispAnalyzerEnergy->initialize( fTLRunParameter->fEnergyReconstruction_BDTFileName, "TMVABDT", "BDTDispEnergy" );
    }


    if( fDebug )
    {
        cout << "VTableLookupDataHandler::setInputFile() END" << endl;
    }

    return fIsMC;
}


/*!

    set data output file and define output tree

        iOutput output file name
        iOption 'RECREATE' or 'UPDATE'

*/
bool VTableLookupDataHandler::setOutputFile( string iOutput, string iOption, string tablefile )
{
    foutputfile = iOutput;

    if( fNTel == 0 )
    {
        cout << "VTableLookupDataHandler::setOutputFile error: no telescopes" << endl;
        exit( EXIT_FAILURE );
    }

    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        if( foutputfile == finputfile[i] && iOption == "recreate" )
        {
            cout << "VTableLookupDataHandler::setOutputFile error: can't overwrite inputfile" << endl;
            cout << "\t" << finputfile[i] << endl;
            exit( EXIT_FAILURE );
        }
    }

    // open output file
    fOutFile = new TFile( foutputfile.c_str(), iOption.c_str() );
    if( fOutFile->IsZombie() )
    {
        cout << "VTableLookupDataHandler::setOutputFile error while opening output file " << foutputfile << "\t" << iOption << endl;
        exit( EXIT_FAILURE );
    }
    // define output tree
    char iTT[2000];

    if( fEventDisplayFileFormat < 6 )
    {
        sprintf( iTT, "MSWC and energy lookup results (%s) VERSION %d", tablefile.c_str(), fEventDisplayFileFormat + 1 );
    }
    else
    {
        sprintf( iTT, "MSWC and energy lookup results (%s) VERSION %d", tablefile.c_str(), fEventDisplayFileFormat );
    }
    fOTree = new TTree( "data", iTT );
    fOTree->SetMaxTreeSize( 1000 * Long64_t( 2000000000 ) );

    fOTree->Branch( "runNumber", &runNumber, "runNumber/I" );
    fOTree->Branch( "eventNumber", &eventNumber, "eventNumber/I" );
    fOTree->Branch( "MJD", &MJD, "MJD/I" );
    fOTree->Branch( "Time",  &time,  "Time/D" );
    sprintf( iTT, "TelElevation[%d]/D", fNTel );
    fOTree->Branch( "TelElevation", fTelElevation, iTT );
    sprintf( iTT, "TelAzimuth[%d]/D", fNTel );
    fOTree->Branch( "TelAzimuth", fTelAzimuth, iTT );
    fOTree->Branch( "ArrayPointing_Elevation", &fArrayPointing_Elevation, "ArrayPointing_Elevation/F" );
    fOTree->Branch( "ArrayPointing_Azimuth", &fArrayPointing_Azimuth, "ArrayPointing_Azimuth/F" );
    fOTree->Branch( "WobbleN", &fWobbleN, "WobbleN/D" );
    fOTree->Branch( "WobbleE", &fWobbleE, "WobbleE/D" );

    fOTree->Branch( "LTrig", &LTrig, "LTrig/l" );
    fOTree->Branch( "NTrig", &fNTrig, "NTrig/i" );
    fOTree->Branch( "NImages", &fNImages, "NImages/I" );
    fOTree->Branch( "ImgSel", &fImgSel, "ImgSel/l" );

    // MC parameters
    if( fIsMC )
    {
        fOTree->Branch( "MCprimary", &fMCPrimary, "MCprimary/I" );
        fOTree->Branch( "MCe0", &fMCEnergy, "MCe0/D" );
        fOTree->Branch( "MCxcore", &fMCxcore, "MCxcore/D" );
        fOTree->Branch( "MCycore", &fMCycore, "MCycore/D" );
        if( !fShortTree )
        {
            fOTree->Branch( "MCxcore_SC", &fMCxcore_SC, "MCxcore_SC/D" );
        }
        if( !fShortTree )
        {
            fOTree->Branch( "MCycore_SC", &fMCycore_SC, "MCycore_SC/D" );
        }
        if( !fShortTree )
        {
            fOTree->Branch( "MCxcos", &fMCxcos, "MCxcos/D" );
        }
        if( !fShortTree )
        {
            fOTree->Branch( "MCycos", &fMCycos, "MCycos/D" );
        }
        fOTree->Branch( "MCaz", &fMCaz, "MCaz/D" );
        fOTree->Branch( "MCze", &fMCze, "MCze/D" );
        fOTree->Branch( "MCxoff", &fMCxoff, "MCxoff/D" );
        fOTree->Branch( "MCyoff", &fMCyoff, "MCyoff/D" );
        fOTree->Branch( "MCCorsikaRunID", &fMCCorsikaRunID, "MCCorsikaRunID/I" );
        fOTree->Branch( "MCCorsikaShowerID", &fMCCorsikaShowerID, "MCCorsikaShowerID/I" );
        fOTree->Branch( "MCFirstInteractionHeight", &fMCFirstInteractionHeight, "MCFirstInteractionHeight/F" );
        fOTree->Branch( "MCFirstInteractionDepth", &fMCFirstInteractionDepth, "MCFirstInteractionDepth/F" );
        fOTree->Branch( "MCR", fR_short_MC, "MCR[NImages]/D" );
    }

    // telescope type related variables
    fOTree->Branch( "ImgSel_list",  fImgSel_list_short, "ImgSel_list[NImages]/i" );
    fOTree->Branch( "NTtype", &fNTelTypes, "NTtype/I" );
    fOTree->Branch( "TtypeID", fTtypeID, "TtypeID[NTtype]/l" );
    fOTree->Branch( "NImages_Ttype", NImages_Ttype, "NImages_Ttype[NTtype]/i" );

    fOTree->Branch( "img2_ang", &fimg2_ang, "img2_ang/D" );
    fOTree->Branch( "RecID", &fMethod, "RecID/I" );
    fOTree->Branch( "Ze", &fZe, "Ze/D" );
    fOTree->Branch( "Az", &fAz, "Az/D" );
    if( !fShortTree )
    {
        fOTree->Branch( "ra", &fRA, "ra/D" );
    }
    if( !fShortTree )
    {
        fOTree->Branch( "dec", &fDec, "dec/D" );
    }
    fOTree->Branch( "Xoff", &fXoff, "Xoff/D" );
    fOTree->Branch( "Yoff", &fYoff, "Yoff/D" );
    fOTree->Branch( "Xoff_derot", &fXoff_derot, "Xoff_derot/D" );
    fOTree->Branch( "Yoff_derot", &fYoff_derot, "Yoff_derot/D" );
    fOTree->Branch( "stdS", &fstdS, "stdS/D" );
    if( !fShortTree )
    {
        fOTree->Branch( "theta2", &ftheta2, "theta2/D" );
    }
    if( !fShortTree )
    {
        fOTree->Branch( "theta2_All", &ftheta2_All, "theta2_All[25]/D" );
    }
    fOTree->Branch( "Xcore", &fXcore, "Xcore/D" );
    fOTree->Branch( "Ycore", &fYcore, "Ycore/D" );
    if( !fShortTree )
    {
        fOTree->Branch( "Xcore_SC", &fXcore_SC, "Xcore_SC/D" );
    }
    if( !fShortTree )
    {
        fOTree->Branch( "Ycore_SC", &fYcore_SC, "Ycore_SC/D" );
    }
    fOTree->Branch( "stdP", &fstdP, "stdP/D" );
    fOTree->Branch( "Chi2", &fchi2, "Chi2/D" );

    fOTree->Branch( "meanPedvar_Image", &fmeanPedvar_Image, "meanPedvar_Image/F" );

    // image parameters (per reconstructed telescope)
    fOTree->Branch( "ntubes", fntubes_short, "ntubes[NImages]/I" );
    fOTree->Branch( "dist", fdist_short, "dist[NImages]/F" );
    fOTree->Branch( "fui", ffui_short, "fui[NImages]/F" );
    fOTree->Branch( "size", fsize_short, "size[NImages]/F" );
    fOTree->Branch( "loss", floss_short, "loss[NImages]/F" );
    fOTree->Branch( "width", fwidth_short, "width[NImages]/F" );
    fOTree->Branch( "dwidth", fdwidth_short, "dwidth[NImages]/F" );
    fOTree->Branch( "length", flength_short, "length[NImages]/F" );
    fOTree->Branch( "dlength", fdlength_short, "dlength[NImages]/F" );
    fOTree->Branch( "cen_x", fcen_x_short, "cen_x[NImages]/F" );
    fOTree->Branch( "dcen_x", fdcen_x_short, "dcen_x[NImages]/F" );
    fOTree->Branch( "cen_y", fcen_y_short, "cen_y[NImages]/F" );
    fOTree->Branch( "dcen_y", fdcen_y_short, "dcen_y[NImages]/F" );
    fOTree->Branch( "cosphi", fcosphi_short, "cosphi[NImages]/F" );
    fOTree->Branch( "sinphi", fsinphi_short, "sinphi[NImages]/F" );
    fOTree->Branch( "dphi", fdphi_short, "dphi[NImages]/F" );
    fOTree->Branch( "tgrad_x", ftgrad_x_short, "tgrad_x[NImages]/F" );
    fOTree->Branch( "asym", fasym_short, "asym[NImages]/F" );
    fOTree->Branch( "Fitstat", fFitstat_short, "Fitstat[NImages]/I" );

    // image parameters (for all telescopes)
    if( !fShortTree )
    {
        sprintf( iTT, "meanPedvar_ImageT[%d]/F", fNTel );
        fOTree->Branch( "meanPedvar_ImageT", fmeanPedvar_ImageT, iTT );
        sprintf( iTT, "size2[%d]/D", fNTel );
        fOTree->Branch( "size2", fsize2, iTT );
        sprintf( iTT, "fracLow[%d]/D", fNTel );
        fOTree->Branch( "fracLow", ffracLow, iTT );
        sprintf( iTT, "max1[%d]/D", fNTel );
        fOTree->Branch( "max1", fmax1, iTT );
        sprintf( iTT, "max2[%d]/D", fNTel );
        fOTree->Branch( "max2", fmax2, iTT );
        sprintf( iTT, "max3[%d]/D", fNTel );
        fOTree->Branch( "max3", fmax3, iTT );
        sprintf( iTT, "maxindex1[%d]/I", fNTel );
        fOTree->Branch( "maxindex1", fmaxindex1, iTT );
        sprintf( iTT, "maxindex2[%d]/I", fNTel );
        fOTree->Branch( "maxindex2", fmaxindex2, iTT );
        sprintf( iTT, "maxindex3[%d]/I", fNTel );
        fOTree->Branch( "maxindex3", fmaxindex3, iTT );
        sprintf( iTT, "nsat[%d]/s", fNTel );
        fOTree->Branch( "nsat", fnsat, iTT );
        sprintf( iTT, "nlowgain[%d]/s", fNTel );
        fOTree->Branch( "nlowgain", fnlowgain, iTT );
        sprintf( iTT, "alpha[%d]/D", fNTel );
        fOTree->Branch( "alpha", falpha, iTT );
        sprintf( iTT, "los[%d]/D", fNTel );
        fOTree->Branch( "los", flos, iTT );
        sprintf( iTT, "cen_x[%d]/D", fNTel );
        fOTree->Branch( "cen_x", fcen_x, iTT );
        sprintf( iTT, "cen_y[%d]/D", fNTel );
        fOTree->Branch( "cen_y", fcen_y, iTT );
        sprintf( iTT, "cosphi[%d]/D", fNTel );
        fOTree->Branch( "cosphi", fcosphi, iTT );
        sprintf( iTT, "sinphi[%d]/D", fNTel );
        fOTree->Branch( "sinphi", fsinphi, iTT );
        sprintf( iTT, "tchisq_x[%d]/D", fNTel );
        fOTree->Branch( "tchisq_x", ftchisq_x, iTT );
    }
    //    if( fTreeWithParameterErrors )
    {
        sprintf( iTT, "dcen_x[%d]/D", fNTel );
        fOTree->Branch( "dcen_x", fdcen_x, iTT );
        sprintf( iTT, "dcen_y[%d]/D", fNTel );
        fOTree->Branch( "dcen_y", fdcen_y, iTT );
        sprintf( iTT, "dlength[%d]/D", fNTel );
        fOTree->Branch( "dlength", fdlength, iTT );
        sprintf( iTT, "dwidth[%d]/D", fNTel );
        fOTree->Branch( "dwidth", fdwidth, iTT );
        sprintf( iTT, "dphi[%d]/D", fNTel );
        fOTree->Branch( "dphi", fdphi, iTT );
    }
    fOTree->Branch( "DispNImages", &fnxyoff, "DispNImages/i" );
    fOTree->Branch( "DispXoff_T", fXoff_T, "DispXoff_T[DispNImages]/F" );
    fOTree->Branch( "DispYoff_T", fYoff_T, "DispYoff_T[DispNImages]/F" );
    fOTree->Branch( "DispWoff_T", fWoff_T, "DispWoff_T[DispNImages]/F" );
    fOTree->Branch( "DispAbsSumWeigth", &fDispAbsSumWeigth, "fDispAbsSumWeigth/F" );
    fOTree->Branch( "Disp_T", fDoff_T, "Disp_T[DispNImages]/F" );
    fOTree->Branch( "DispTelList_T", fToff_T, "DispTelList_T[DispNImages]/i" );
    fOTree->Branch( "DispDiff", &fDispDiff, "DispDiff/D" );
    fOTree->Branch( "Xoff_intersect", &fXoff_intersect, "Xoff_intersect/F" );
    fOTree->Branch( "Yoff_intersect", &fYoff_intersect, "Yoff_intersect/F" );
    fOTree->Branch( "cross", fcross_short, "cross[NImages]/F" );
    fOTree->Branch( "crossO", fcrossO_short, "crossO[NImages]/F" );
    if( fIsDeepLearner )
    {
        fOTree->Branch( "dl_gammaness", &dl_gammaness, "dl_gammaness/D" );
        fOTree->Branch( "dl_isGamma", &dl_isGamma, "dl_isGamma/O" );
    }

    fOTree->Branch( "R", fR_short, "R[NImages]/D" );
    if( !fShortTree )
    {
        sprintf( iTT, "MSCWT[%d]/D", fNTel );
        fOTree->Branch( "MSCWT", ftmscw, iTT );
        sprintf( iTT, "MSCLT[%d]/D", fNTel );
        fOTree->Branch( "MSCLT", ftmscl, iTT );
        if( fTLRunParameter && fTLRunParameter->fUsetimeGradientLookupTables )
        {
            sprintf( iTT, "MSCTT[%d]/D", fNTel );
            fOTree->Branch( "MSCTT", ftmsct, iTT );
        }
        sprintf( iTT, "MSCWTSigma[%d]/F", fNTel );
        fOTree->Branch( "MSCWTSigma", ftmscw_sigma, iTT );
        sprintf( iTT, "MSCLTSigma[%d]/F", fNTel );
        fOTree->Branch( "MSCLTSigma", ftmscl_sigma, iTT );
        if( fTLRunParameter && fTLRunParameter->fUsetimeGradientLookupTables )
        {
            sprintf( iTT, "MSCTTSigma[%d]/F", fNTel );
            fOTree->Branch( "MSCTTSigma", ftmsct_sigma, iTT );
        }
    }
    sprintf( iTT, "NMSCW/I" );
    fOTree->Branch( "NMSCW", &fnmscw, iTT );
    fOTree->Branch( "ES", fES_short, "ES[NImages]/D" );

    sprintf( iTT, "MSCW/D" );
    fOTree->Branch( "MSCW", &fmscw, iTT );
    sprintf( iTT, "MSCL/D" );
    fOTree->Branch( "MSCL", &fmscl, iTT );
    sprintf( iTT, "MWR/F" );
    fOTree->Branch( "MWR", &fmwr, iTT );
    sprintf( iTT, "MLR/F" );
    fOTree->Branch( "MLR", &fmlr, iTT );
    if( fTLRunParameter && fTLRunParameter->fUsetimeGradientLookupTables )
    {
        sprintf( iTT, "MSCT/D" );
        fOTree->Branch( "MSCT", &fmsct, iTT );
    }
    sprintf( iTT, "ErecS/D" );
    fOTree->Branch( "ErecS", &fenergyS, iTT );
    sprintf( iTT, "EChi2S/D" );
    fOTree->Branch( "EChi2S", &fechi2S, iTT );
    sprintf( iTT, "dES/D" );
    fOTree->Branch( "dES", &fdES, iTT );
    sprintf( iTT, "dESabs/D" );
    fOTree->Branch( "dESabs", &feAbsError, iTT );
    sprintf( iTT, "NErecST/I" );
    fOTree->Branch( "NErecST", &fnenergyT, iTT );
    fOTree->Branch( "ErecQL", &fenergyQL, "ErecQL/I" );

    sprintf( iTT, "EmissionHeight/F" );
    fOTree->Branch( "EmissionHeight", &fEmissionHeightMean, iTT );
    sprintf( iTT, "EmissionHeightChi2/F" );
    fOTree->Branch( "EmissionHeightChi2", &fEmissionHeightChi2, iTT );
    sprintf( iTT, "NTelPairs/i" );
    fOTree->Branch( "NTelPairs", &fNTelPairs, iTT );
    fOTree->Branch( "SizeSecondMax", &fSizeSecondMax, "SizeSecondMax/D" );

    sprintf( iTT, "fEmissionHeightT[NTelPairs]/F" );
    if( !fShortTree )
    {
        fOTree->Branch( "EmissionHeightT", fEmissionHeightT, iTT );
    }
    for( unsigned int i = 0; i < getMaxNbrTel(); i++ )
    {
        fEmissionHeightT[i] = -99.;
    }
    if( fTLRunParameter->fWritePixelLists )
    {
        fOTree->Branch( "PixelListN", PixelListN, "PixelListN[NImages]/i" );
        fOTree->Branch( "PixelListNPixelNN", &PixelListNPixelNN, "PixelListNPixelNN/i" );
        fOTree->Branch( "PixelID", PixelID, "PixelID[PixelListNPixelNN]/i" );
        fOTree->Branch( "PixelType", PixelType, "PixelType[PixelListNPixelNN]/i" );
        fOTree->Branch( "PixelIntensity", PixelIntensity, "PixelIntensity[PixelListNPixelNN]/F" );
        fOTree->Branch( "PixelTimingT0", PixelTimingT0, "PixelTimingT0[PixelListNPixelNN]/F" );
        fOTree->Branch( "PixelPE", PixelPE, "PixelPE[PixelListNPixelNN]/F" );
    }

    readRunParameter();

    return true;
}


/*
 * read and update run parameters from eventdisplay file
 *
 * Note: read all run parameter from first non-Zombie file
 *
 */
bool VTableLookupDataHandler::readRunParameter()
{
    if( fEventDisplayFileFormat > 1 )
    {
        // get list of files in chain
        TObjArray* fileElements = fTshowerpars->GetListOfFiles();
        TChainElement* chEl = 0;
        TIter next( fileElements );
        chEl = ( TChainElement* )next();

        TFile ifInput( chEl->GetTitle() );
        if( !ifInput.IsZombie() )
        {
            cout << "reading eventdisplay run parameters from " << ifInput.GetName() << endl;
            TNamed* iR = ( TNamed* )ifInput.Get( "runparameter" );
            if( iR && fOutFile )
            {
                fOutFile->cd();
                iR->Write();
            }
            VEvndispRunParameter* iPar = ( VEvndispRunParameter* ) ifInput.Get( "runparameterV2" );
            VEvndispReconstructionParameter* iERecPar = ( VEvndispReconstructionParameter* )ifInput.Get( "EvndispReconstructionParameter" );
            VMonteCarloRunHeader* iMC = ( VMonteCarloRunHeader* )ifInput.Get( "MC_runheader" );
            if( iMC )
            {
                fTLRunParameter->ze = iMC->getMeanZenithAngle_Deg();
            }
            if( iPar )
            {
                if( fTLRunParameter->fTelToAnalyse.size() > 0 )
                {
                    iPar->fTelToAnalyze = fTLRunParameter->fTelToAnalyse;
                }
                else if( iERecPar )
                {
                    //copied from VTableLookupRunParameter.cpp
                    vector< unsigned int > iRunParT = iPar->fTelToAnalyze;
                    vector< unsigned int > iTelToAnalyze;
                    // this works only if number of telescopes = number of telescope types
                    // (e.g. the VERITAS case)
                    if( iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )
                            && iPar->fNTelescopes == iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )->fLocalUseImage.size() )
                    {
                        for( unsigned int i = 0; i < iRunParT.size(); i++ )
                        {
                            if( iRunParT[i] < iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )->fLocalUseImage.size()
                                    && iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )->fLocalUseImage[iRunParT[i]] )
                            {
                                iTelToAnalyze.push_back( iRunParT[i] );
                            }
                        }
                    }
                    iPar->fTelToAnalyze = iTelToAnalyze;
                }
                if( fOutFile )
                {
                    fOutFile->cd();
                    // update instrument epoch in evendisp run parameters
                    // (might have been changed since the evndisp analysis)
                    if( fTLRunParameter->fUpdateInstrumentEpoch )
                    {
                        cout << "Evaluating instrument epoch (";
                        cout << "was: " << iPar->getInstrumentEpoch( false );
                        cout << ", is: " << iPar->getInstrumentEpoch( false, true );
                        cout << ")" << endl;
                        cout << "Evaluating atmosphere ID (";
                        cout << "was: " << iPar->getAtmosphereID( false );
                        cout << ", is: " << iPar->getAtmosphereID( true );
                        cout << ")" << endl;
                    }
                    iPar->Write();
                }
            }
            ///////////////////////////
            // read parameters needed for (simple) stereo reconstruction
            if( iERecPar && iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method ) )
            {
                // minimum angle set as command line parameter
                if( fTLRunParameter && fTLRunParameter->fRerunStereoReconstruction_minAngle > 0. )
                {
                    fSSR_AxesAngles_min = fTLRunParameter->fRerunStereoReconstruction_minAngle;
                }
                // use minimum angle from evndisp analysis
                else
                {
                    fSSR_AxesAngles_min = iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )->fAxesAngles_min;
                }
                fSSR_NImages_min    = iERecPar->getReconstructionParameterData( fTLRunParameter->rec_method )->fNImages_min;
                if( fTLRunParameter->fRerunStereoReconstruction )
                {
                    cout << "\t quality cuts in stereo reconstruction: ";
                    cout << "number of images >= " << fSSR_NImages_min;
                    cout << ", angdiff > " << fSSR_AxesAngles_min << " deg" << endl;
                }
            }
            ///////////////////////////
            // simulated energy spectrum
            // and zenith angle
            if( iMC )
            {
                setMCMinEnergy( iMC->E_range[0] );
                setMCMaxEnergy( iMC->E_range[1] );
                setMCSpectralIndex( -1 * iMC->spectral_index );
                fTLRunParameter->ze = iMC->getMeanZenithAngle_Deg();
            }
            ifInput.Close();
            if( fOutFile )
            {
                fOutFile->cd();
            }
        }
    }

    return true;
}

void VTableLookupDataHandler::printCutStatistics()
{
    cout << "---------------------------------------------------------------------------------------------------" << endl;
    cout << "Cut statistics: " << endl;
    if( fNStats_All == 0 )
    {
        cout << "\t no events..." << endl;
        return;
    }
    unsigned int nTOT = fNStats_All;

    cout << "\t number of events considered: \t\t" << fNStats_All << endl;
    nTOT -= fNStats_NImagesCut;
    cout << "\t removed by >= " << fTLRunParameter->fTableFillingCut_NImages_min  << " images: \t\t\t" << fNStats_NImagesCut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_NImagesCut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    nTOT = nTOT - fNStats_Chi2Cut + fNStats_NImagesCut;
    cout << "\t removed by Chi2 >=0:   \t\t\t" << fNStats_Chi2Cut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_Chi2Cut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    cout << "\t number of reconstructed events:   \t" << fNStats_Rec;
    cout << " (fraction of reconstructed events: " << ( float )fNStats_Rec / ( float )fNStats_All << "; " << nTOT << ")" << endl;

    nTOT -= fNStats_CoreErrorCut;
    cout << "\t removed by cut on core misreconstruction: \t\t" << fNStats_CoreErrorCut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_CoreErrorCut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    nTOT = nTOT - fNStats_WobbleCut + fNStats_CoreErrorCut;
    cout << "\t removed by wobble cut (<" << fTLRunParameter->fTableFillingCut_WobbleCut_max << "): \t\t\t" << fNStats_WobbleCut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_WobbleCut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    nTOT = nTOT - fNStats_WobbleMinCut + fNStats_WobbleCut;
    cout << "\t removed by MC wobble min cut (>" << fMC_distance_to_cameracenter_min << "): \t\t" << fNStats_WobbleMinCut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_WobbleMinCut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    nTOT -= fNStats_WobbleMaxCut;
    cout << "\t removed by wobble max cut (<" << fMC_distance_to_cameracenter_max << "): \t\t" << fNStats_WobbleMaxCut;
    cout << " (fraction removed/# of events left: " << ( float )fNStats_WobbleMaxCut / ( float )fNStats_All << "; " << nTOT << ")" << endl;
    cout << "---------------------------------------------------------------------------------------------------" << endl;
}

/*!
  write everything to disk
*/
bool VTableLookupDataHandler::terminate( TNamed* iM )
{
    printCutStatistics();

    if( fOutFile )
    {
        cout << "writing data to " << fOutFile->GetName() << endl;
        fOutFile->cd();

        cout << endl << "\t total number of events in output tree: " << fOTree->GetEntries() << endl << endl;
        fOTree->Write( "", TObject::kOverwrite );

        if( iM )
        {
            cout << "\t writing table lookup run parameter" << endl;
            iM->Write();
        }
        else
        {
            cout << "\t no table lookup run parameter to write" << endl;
        }

        if( fTtelconfig )
        {
            cout << "\t writing telescope configuration" << endl;
            copy_telconfig();
        }
        else
        {
            cout << "\t no telescope configuration to write" << endl;
        }
        fOutFile->cd();

        if( fIsMC )
        {
            cout << "\t writing MC debug histograms" << endl;
            hisList->Write();
        }
        // see if there is a dead time object on file, if not: write the one filled here
        // (note: at this stage, the scalars cannot be used and the dead time might be
        //         underestimated)
        if( !fIsMC )
        {
            writeDeadTimeHistograms();
        }

        // copy MC tree
        // (not default, as this is a large tree with
        // 1 entry per simulated event)
        if( fIsMC )
        {
            bool iMCTree_exists = copyMCRunheader();
            if( bWriteMCPars && iMCTree_exists )
            {
                copyTree_from_evndispFile();
            }
            copyMCHistograms();
        }
        // try and copy a deep learning tree
        // (usually not there)
        // copyTree_from_evndispFile( "data_DL" );


        ///////////////////////////////////////////////////////////////////////////
        // copy TTree 'pointingDataReduced' and 'deadPixelRegistry' from evndisp.<>.root to mscw.<>.root
        if( finputfile.size() > 1 && !fIsMC )
        {
            cout << "Warning, VTableLookupDataHandler->finputfile.size() isn't 1, not sure which input file to copy TTree 'pointingDataReduced' from, copying from file finputfile[0]:" << finputfile[0] << endl;
        }
        // not sure why we don't want to do this for MC
        if( finputfile.size() > 0 && !fIsMC )
        {
            TFile* inpMscwFile = new TFile( finputfile[0].c_str(), "READ" ) ;
            fOutFile->cd();
            TTree* iTree       = ( TTree* )inpMscwFile->Get( "pointingDataReduced" );
            if( iTree )
            {
                TTree* newtree     = iTree->CloneTree();
                if( newtree )
                {
                    newtree->Write();
                }
                else
                {
                    cout << "VTableLookupDataHandler::terminate Warning: Unable to clone tree " << iTree->GetName() << endl;
                }
            }
            else
            {
                cout << "VTableLookupDataHandler::terminate Warning: Unable to find tree pointingDataReduced in file " << inpMscwFile->GetName() << endl;
            }

            TTree* jTree = ( TTree* )inpMscwFile->Get( "deadPixelRegistry" ) ;
            // deadPixelRegistry may not exist, only try to copy it if it's there
            if( jTree )
            {
                fOutFile->cd() ;
                TTree* newtree2 = jTree->CloneTree() ;
                newtree2->Write() ;
            }

        }
        else if( !fIsMC )
        {
            cout << "Warning, VTableLookupDataHandler->finputfile has size 0, unable to copy TTree 'pointingDataReduced' to file " << fOutFile->GetName() << endl;
        }


        fOutFile->Close();
        cout << "...outputfile closed" << endl;
        cout << "(" << fOutFile->GetName() << ")" << endl;
    }

    return true;
}

void VTableLookupDataHandler::writeDeadTimeHistograms()
{
    if( finputfile.size() > 1 )
    {
        cout << "VTableLookupDataHandler::writeDeadTimeHistograms() error: ";
        cout << "analysis of several files at once not allowed " << endl;
        cout << "(dead times will be wrong)" << endl;
        exit( EXIT_FAILURE );
    }

    // use chain to get list of files
    TChain iTel( "telconfig" );
    int iNFil = 0;
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        iNFil += iTel.Add( finputfile[i].c_str() );
    }

    if( iNFil > 0 )
    {
        TFile* f = iTel.GetFile();
        if( f )
        {
            TDirectoryFile* iDeadtimeDirectory = ( TDirectoryFile* )f->Get( "deadTimeHistograms" );
            if( iDeadtimeDirectory )
            {
                VDeadTime* iDeadTime = new VDeadTime();
                iDeadTime->readHistograms( iDeadtimeDirectory );
                iDeadTime->calculateDeadTime();
                iDeadTime->printDeadTime();
                iDeadTime->writeHistograms();
            }
        }
    }
}

/*
 * copy tel_config tree from first evndisp file into
 * mscw output file
 *
 * copy only telescopes which are selected
 * with the command line parameter
 *
 */
void VTableLookupDataHandler::copy_telconfig()
{
    TChain iMC( "telconfig" );
    int iNFil = 0;
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        iNFil += iMC.Add( finputfile[i].c_str() );
    }

    if( iNFil > 0 )
    {
        TFile* f = iMC.GetFile();
        if( f && fTLRunParameter )
        {
            if( f->Get( "telconfig" ) )
            {
                TTree* t = ( TTree* )f->Get( "telconfig" );
                Ctelconfig* tc = new Ctelconfig( t );
                fOutFile->cd();
                TTree* n = tc->fChain->CloneTree( 0 );
                if( fTLRunParameter && tc->fChain->GetEntries() != ( int )fTLRunParameter->fTelToAnalyzeData.size() )
                {
                    cout << "VTableLookupDataHandler::copy_telconfig() error: ";
                    cout << "mismatch between telecope vector and telconfig tree";
                    return;
                }
                // get total number of valid telescopes;
                UInt_t i_ntel = 0;
                for( unsigned int i = 0; i < fTLRunParameter->fTelToAnalyzeData.size(); i++ )
                {
                    if( fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
                    {
                        i_ntel++;
                    }
                }

                // copy only selected telescopes
                for( unsigned int i = 0; i < fTLRunParameter->fTelToAnalyzeData.size(); i++ )
                {
                    tc->GetEntry( i );
                    if( fTLRunParameter->fTelToAnalyzeData[i] && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
                    {
                        // fix total number of telescopes
                        tc->NTel = i_ntel;

                        n->Fill();
                    }
                }
                n->Write();
            }
        }
    }
}

/*
   copy MC run header from first file

   OBS: assume that all run headers in all files are the same!!!!

   check additionally if MCpars tree exists
*/
bool VTableLookupDataHandler::copyMCRunheader()
{
    // use chain to get list of files
    TChain iTel( "telconfig" );
    int iNFil = 0;
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        iNFil += iTel.Add( finputfile[i].c_str() );
    }

    if( iNFil > 0 )
    {
        TFile* f = iTel.GetFile();
        if( f )
        {
            if( ( VMonteCarloRunHeader* )f->Get( "MC_runheader" ) )
            {
                fOutFile->cd();
                f->Get( "MC_runheader" )->Write();
                cout << "\t MC run header found and copied" << endl;
            }
            if( f->Get( "MCpars" ) )
            {
                return true;
            }
        }
    }
    return false;
}

/*
 * copy a tree from the eventdisplay files to the mscw_energy
 * output file
 *
 */
void VTableLookupDataHandler::copyTree_from_evndispFile( string iTreeName )
{
    TChain iMC( iTreeName.c_str() );
    int iNFil = 0;
    for( unsigned int i = 0; i < finputfile.size(); i++ )
    {
        iNFil += iMC.Add( finputfile[i].c_str() );
    }

    if( iNFil > 0 && iMC.GetEntries() > 0 )
    {
        cout << "\t copying " << iTreeName;
        cout << " with " << iMC.GetEntries() << " entries..." << flush;
        iMC.Merge( fOutFile, 0, "keep" );
        cout << "done " << endl;
    }
}

void VTableLookupDataHandler::copyMCHistograms()
{
    vector< double > i_az_min;
    vector< double > i_az_max;
    vector< double > i_spectral_index = fTLRunParameter->fAddMC_spectral_index;;
    if( i_spectral_index.size() == 1 )
    {
        cout << "\t VTableLookupDataHandler::copyMCHistograms: reducing spectral index range to one value: ";
        cout << i_spectral_index[0] << endl;
    }
    VEffectiveAreaCalculatorMCHistograms* iMC_his = 0;
    if( fTshowerpars )
    {
        // loop over all files in chain (might be many) and add up MC histograms
        // (histograms are needed for effective area calculation)
        TObjArray* fileElements = fTshowerpars->GetListOfFiles();
        if( !fileElements )
        {
            cout << "VTableLookupDataHandler::copyMCHistograms(): no list of files found" << endl;
            return;
        }
        TChainElement* chEl = 0;
        TIter next( fileElements );
        unsigned int z = 0;
        unsigned int i_failure = 0;
        unsigned int i_all = 0;
        while( ( chEl = ( TChainElement* )next() ) )
        {
            TFile* ifInput = new TFile( chEl->GetTitle() );
            if( !ifInput->IsZombie() )
            {
                // first file: all following MC histograms are added to this one
                // (file is not closed)
                if( z == 0 )
                {
                    iMC_his = ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" );
                    if( iMC_his && i_spectral_index.size() > 0 )
                    {
                        iMC_his->matchDataVectors( i_az_min, i_az_max, i_spectral_index );
                    }
                }
                else
                {
                    if( iMC_his )
                    {
                        VEffectiveAreaCalculatorMCHistograms* iMC_his_temp =
                            ( VEffectiveAreaCalculatorMCHistograms* )ifInput->Get( "MChistos" );
                        if( iMC_his_temp && i_spectral_index.size() > 0 )
                        {
                            iMC_his_temp->matchDataVectors( i_az_min, i_az_max, i_spectral_index );
                        }
                        if( iMC_his_temp )
                        {
                            iMC_his->add( iMC_his_temp );
                        }
                    }
                    ifInput->Close();
                }
                z++;
            }
            else
            {
                cout << "VTableLookupDataHandler::copyMCHistograms(): warning: no file found for MC histograms" << endl;
                i_failure++;
            }
            i_all++;
        }
        cout << "\t Search " << i_all << " files for MC histograms (found: " << z;
        cout << ", failed: " << i_failure << ")" << endl;
        if( iMC_his && fOutFile )
        {
            cout << "\t writing MC histograms" << endl;
            iMC_his->print();
            fOutFile->cd();
            iMC_his->Write();
        }
    }
}


void VTableLookupDataHandler::reset()
{
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fR[i] = -99.;
        fRTel[i] = -99.;
        fR_short[i] = -99.;
        fcross_short[i] = -99.;
        fcrossO_short[i] = -99.;
        fntubes_short[i] = -99;
        fdist_short[i] = -99.;
        ffui_short[i] = -99.;
        fsize_short[i] = -99.;
        floss_short[i] = -99.;
        fwidth_short[i] = -99.;
        fdwidth_short[i] = -99.;
        flength_short[i] = -99.;
        fdlength_short[i] = -99.;
        fcen_x_short[i] = -99.;
        fcen_y_short[i] = -99.;
        fdcen_x_short[i] = -99.;
        fdcen_y_short[i] = -99.;
        fcosphi_short[i] = -99.;
        fsinphi_short[i] = -99.;
        fdphi_short[i] = -99.;
        ftgrad_x_short[i] = -99.;
        fasym_short[i] = -99.;
        fFitstat_short[i] = -99;
        fR_MC[i] = -99.;
        fR_short_MC[i] = -99.;
        fES[i] = -99.;
        fES_short[i] = -99.;
        ftmscl[i] = -99.;
        ftmscw[i] = -99.;
        ftmsct[i] = -99.;
        ftmscw_sigma[i] = -99.;
        ftmscl_sigma[i] = -99.;
        ftmsct_sigma[i] = -99.;
        fXoff_T[i] = -99.;
        fYoff_T[i] = -99.;
        fWoff_T[i] = -99.;
        fDoff_T[i] = -99.;
        fToff_T[i] = 0;
        if( fTLRunParameter->fWritePixelLists )
        {
            PixelListN[i] = 0;
        }

    }
    PixelListNPixelNN = 0;
    fnxyoff = 0;
    fnmscw = 0;
    fnenergyT = 0;
    fenergyQL = -1;
    fmscl = -99.;
    fmscw = -99.;
    fmwr  = -99.;
    fmlr  = -99.;
    fenergyS = -99.;
    fechi2S = 1.e-10;
    feAbsError = -99.;
    fdES = -99.;
    fmsct = -99.;
    fXoff = -99.;
    fYoff = -99.;
    fXoff_derot = -99.;
    fYoff_derot = -99.;
    fXoff_intersect = -99.;
    fYoff_intersect = -99.;
    fXoff_edisp = -99.;
    fYoff_edisp = -99.;
    fXcore = -99.;
    fYcore = -99.;
    fstdP = -99.;

    fEmissionHeightMean = -99.;
    fEmissionHeightChi2 = -99.;

    fmeanPedvar_Image = 0.;

    resetImageParameters();
}


/*!
  calculate distances between telescopes and reconstructed shower core

*/
void VTableLookupDataHandler::calcDistances()
{
    // check for successful reconstruction
    for( unsigned int tel = 0; tel < fNTel; tel++ )
    {
        fR_MC[tel] = VUtilities::line_point_distance( fMCycore, -1.*fMCxcore, 0., fMCze, fMCaz, fTelY[tel], -1.*fTelX[tel], fTelZ[tel] );
        if( fImgSel_list[tel] && fZe >= 0. && fXcore > -9998. && fYcore > -9998. )
        {
            fR[tel]    = VUtilities::line_point_distance( fYcore, -1.*fXcore, 0., fZe, fAz, fTelY[tel], -1.*fTelX[tel], fTelZ[tel] );
            fRTel[tel] = VUtilities::line_point_distance( fYcore, -1.*fXcore, 0., 90. - fArrayPointing_Elevation, fArrayPointing_Azimuth, fTelY[tel], -1.*fTelX[tel], fTelZ[tel] );
        }
        else
        {
            fR[tel] = -99.;
            fRTel[tel] = -99.;
        }
    }
}



void VTableLookupDataHandler::resetImageParameters()
{
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        resetImageParameters( i );
    }
}


void VTableLookupDataHandler::resetImageParameters( unsigned int i )
{
    if( i > getMaxNbrTel() )
    {
        return;
    }

    fdist[i] = 0.;
    ffui[i] = 0.;
    fsize[i] = 0.;
    fsize2[i] = 0.;
    fweight[i] = 1.;
    fweightDispBDTs[i] = 1.;
    floss[i] = 0.;
    ffracLow[i] = 0.;
    fwidth[i] = 0.;
    flength[i] = 0.;
    fmeanPedvar_ImageT[i] = 0.;
    if( fwrite )
    {
        return;
    }

    fntubes[i] = 0;
    fnsat[i] = 0;
    fnlowgain[i] = 0;
    falpha[i] = 0.;
    flos[i] = 0.;
    fasym[i] = 0.;
    fcen_x[i] = 0.;
    fcen_y[i] = 0.;
    fdcen_x[i] = 0.;
    fdcen_y[i] = 0.;
    fdwidth[i] = 0.;
    fdlength[i] = 0.;
    fdphi[i] = 0.;
    fcosphi[i] = 0.;
    fsinphi[i] = 0.;
    fmax1[i] = 0.;
    fmax2[i] = 0.;
    fmax3[i] = 0.;
    fmaxindex1[i] = 0;
    fmaxindex2[i] = 0;
    fmaxindex3[i] = 0;
    ftgrad_x[i] = 0.;
    ftchisq_x[i] = 0.;
    fFitstat[i] = 0;
    if( fTLRunParameter->fWritePixelLists )
    {
        PixelListN[i] = 0;
    }
}

/*
 *
 * quick test if an event has been successfully
 * reconstructed in eventdisplay
 *
 */
bool VTableLookupDataHandler::isReconstructed( bool iEventCounters )
{
    // use MC parameters for table filling
    if( fTLRunParameter->fTableFilling_useStereoMCParameter
            && fNImages >= ( int )fTLRunParameter->fTableFillingCut_NImages_min )
    {
        return true;
    }
    // require successful reconstruction
    if( fchi2 < 0 && fNImages < ( int )fTLRunParameter->fTableFillingCut_NImages_min )
    {
        if( iEventCounters )
        {
            fNStats_Chi2Cut++;
        }
        return false;
    }
    // require stereo events with the given multiplicity
    // (note that for subarray selection in mscw_energy,
    //  this is only a low limit cut. Further telescopes
    //  can be removed)
    if( fNImages < ( int )fTLRunParameter->fTableFillingCut_NImages_min )
    {
        if( iEventCounters )
        {
            fNStats_NImagesCut++;
        }
        return false;
    }

    return true;
}

/*
 *  calculate emission height
 *
 *  - mean value
 *  - value for each pair
 *  - chi2 for paris
 */
void VTableLookupDataHandler::calcEmissionHeights()
{
    fEmissionHeightCalculator->getEmissionHeight( fcen_x, fcen_y, fsize,
            fArrayPointing_Azimuth, fArrayPointing_Elevation );
    fNTelPairs = fEmissionHeightCalculator->getNTelPairs();
    fEmissionHeightMean = ( float )fEmissionHeightCalculator->getMeanEmissionHeight();
    fEmissionHeightChi2 = ( float )fEmissionHeightCalculator->getMeanEmissionHeightChi2();
    if( fEmissionHeightChi2 <= 0. )
    {
        fEmissionHeightChi2 = 1.e-10;
    }
    for( unsigned int i = 0; i < fNTelPairs; i++ )
    {
        if( i >= getMaxNbrTel() || i >= fEmissionHeightCalculator->getEmissionHeights().size() )
        {
            break;
        }
        fEmissionHeightT[i] = ( float )fEmissionHeightCalculator->getEmissionHeights()[i];
    }
}


void VTableLookupDataHandler::setEventWeightfromMCSpectrum()
{
    fEventWeight = 1.;

    if( TMath::Abs( fSpectralIndex - fMCSpectralIndex ) < 1.e-2 )
    {
        fEventWeight = 1.;
    }
    else if( fSpectralIndex > fMCSpectralIndex )
    {
        double alpha = TMath::Power( fMCMinEnergy, -1.*fMCSpectralIndex ) / TMath::Power( fMCMinEnergy, -1.*fSpectralIndex );

        fEventWeight = alpha * TMath::Power( fMCEnergy, -1.*fSpectralIndex ) / TMath::Power( fMCEnergy, -1.*fMCSpectralIndex );
    }
    else if( fSpectralIndex < fMCSpectralIndex )
    {
        double alpha = TMath::Power( fMCMaxEnergy, -1.*fMCSpectralIndex ) / TMath::Power( fMCMaxEnergy, -1.*fSpectralIndex );

        fEventWeight = alpha * TMath::Power( fMCEnergy, -1.*fSpectralIndex ) / TMath::Power( fMCEnergy, -1.*fMCSpectralIndex );
    }
}


double VTableLookupDataHandler::getZe()
{
    return fZe;
}


/*!

  return copy of fMCEnergy as an array

  (note that mscw modifies these values)

*/
double* VTableLookupDataHandler::getMCEnergyArray()
{
    for( unsigned int i = 0; i < fNTel; i++ )
    {
        fMCEnergyArray[i] = getMCEnergy();
    }
    return fMCEnergyArray;
}


void VTableLookupDataHandler::resetAll()
{
    fEventStatus = true;
    runNumber = 0;
    eventNumber = 0;
    MJD = 0;
    time = 0;
    for( unsigned int i = 0; i < getMaxNbrTel(); i++ )
    {
        fTelElevation[i] = 0.;
        fTelAzimuth[i] = 0.;
    }
    fArrayPointing_Elevation = 0.;
    fArrayPointing_Azimuth = 0.;
    fWobbleN = 0.;
    fWobbleE = 0.;
    fMCPrimary = -99;
    fMCEnergy = 0.;
    fMCxcore = 0.;
    fMCycore = 0.;
    fMCxcore_SC = 0.;
    fMCycore_SC = 0.;
    fMCxcos = 0.;
    fMCycos = 0.;
    fMCaz = 0.;
    fMCze = 0.;
    fMCxoff = 0.;
    fMCyoff = 0.;
    fMCCorsikaRunID = 0;
    fMCCorsikaShowerID = 0;
    fMCFirstInteractionHeight = 0;
    fMCFirstInteractionDepth = 0;
    LTrig = 0;
    fNTrig = 0;
    fNImages = 0;
    fImgSel = 0;
    fNTelTypes = 0;
    for( unsigned int i = 0; i < getMaxNbrTel(); i++ )
    {
        fImgSel_list[i] = false;
        fImgSel_list_short[i] = 0;
        NImages_Ttype[i] = 0;
        fTtypeID[i] = 0;
        fNImages_QCTree[i] = 0;
        fchi2_QCTree[i] = 0.;
    }
    fimg2_ang = 0.;
    fZe = 0.;
    fAz = 0.;
    fRA = 0.;
    fDec = 0.;
    fXoff = 0.;
    fYoff = 0.;
    fXoff_derot = 0.;
    fYoff_derot = 0.;
    fXoff_intersect = 0.;
    fYoff_intersect = 0.;
    fXoff_edisp = 0.;
    fYoff_edisp = 0.;
    fstdS = 0.;
    ftheta2 = 0.;
    fXcore = 0.;
    fYcore = 0.;
    fXcore_SC = 0.;
    fYcore_SC = 0.;
    fstdP = 0.;
    fchi2 = 0.;
    fmeanPedvar_Image = 0.;
    for( unsigned int i = 0; i < getMaxNbrTel(); i++ )
    {
        fdist[i] = 0.;
        ffui[i] = 0.;
        fsize[i] = 0.;
        fsize2[i] = 0.;
        fweight[i] = 1.;
        fweightDispBDTs[i] = 1.;
        fsizeCorr[i] = 0.;
        fsize_telType[i] = 0.;
        floss[i] = 0.;
        ffracLow[i] = 0.;
        fmax1[i] = 0.;
        fmax2[i] = 0.;
        fmax3[i] = 0.;
        fmaxindex1[i] = 0;
        fmaxindex2[i] = 0;
        fmaxindex3[i] = 0;
        fwidth[i] = 0.;
        flength[i] = 0.;
        fntubes[i] = 0;
        fmeanPedvar_ImageT[i] = 0.;
        fnsat[i] = 0;
        fnlowgain[i] = 0;
        falpha[i] = 0.;
        flos[i] = 0.;
        fasym[i] = 0.;
        fcen_x[i] = 0.;
        fcen_y[i] = 0.;
        fcosphi[i] = 0.;
        fsinphi[i] = 0.;
        ftgrad_x[i] = 0.;
        fFitstat[i] = 0;
        ftchisq_x[i] = 0.;
        fR[i] = 0.;
        fRTel[i] = 0.;
        fR_short[i] = 0.;
        fcross_short[i] = 0.;
        fcrossO_short[i] = 0.;
        fntubes_short[i] = -99;
        fdist_short[i] = -99.;
        ffui_short[i] = -99.;
        fsize_short[i] = -99.;
        floss_short[i] = -99.;
        fwidth_short[i] = -99.;
        fdwidth_short[i] = -99.;
        flength_short[i] = -99.;
        fdlength_short[i] = -99.;
        fcen_x_short[i] = -99.;
        fcen_y_short[i] = -99.;
        fdcen_x_short[i] = -99.;
        fdcen_y_short[i] = -99.;
        fcosphi_short[i] = -99.;
        fsinphi_short[i] = -99.;
        fdphi_short[i] = -99.;
        ftgrad_x_short[i] = -99.;
        fasym_short[i] = -99.;
        fFitstat_short[i] = -99;
        fR_telType[i] = 0.;
        fLoss_telType[i] = 0.;
        fDistance_telType[i] = 0.;
        fR_MC[i] = 0.;
        fR_short_MC[i] = 0.;
        fR_telType_MC[i] = 0.;
        fXoff_T[i] = 0.;
        fYoff_T[i] = 0.;
        fWoff_T[i] = 0.;
        fDoff_T[i] = 0.;
        fToff_T[i] = 0;
        ftmscw[i] = 0.;
        ftmscl[i] = 0.;
        ftmsct[i] = 0.;
        ftmscw_sigma[i] = 0.;
        ftmscl_sigma[i] = 0.;
        ftmsct_sigma[i] = 0.;
        fES[i] = 0.;
        fES_short[i] = 0.;
    }
    for( unsigned int i = 0; i < 25; i++ )
    {
        ftheta2_All[i] = 99.;
    }
    fnxyoff = 0;
    fnmscw = 0;
    fnenergyT = 0;
    fenergyQL = -1;
    fmscw = 0.;
    fmscl = 0.;
    fmsct = 0.;
    fmwr  = 0.;
    fmlr  = 0.;
    fenergyS = 0.;
    fechi2S = 1.e-10;
    feAbsError = 0.;
    fdES = 0.;
    fEmissionHeightMean = 0.;
    fEmissionHeightChi2 = 0.;
    fNTelPairs = 0;
    for( unsigned int i = 0; i < getMaxNbrTel(); i++ )
    {
        fEmissionHeightT[i] = 0.;
    }
    fSizeSecondMax = 0.;

    fTotalTime = 0.;
    fTotalTime0 = 0.;

    fMC_distance_to_cameracenter_min = 0.;
    fMC_distance_to_cameracenter_max = 1.e10;
    fDispDiff = 0;
    fDispAbsSumWeigth = 0.;
    // deep learner parameters
    dl_gammaness = 0.;
    dl_isGamma = false;
    fXoff_intersect = 0.;
    fYoff_intersect = 0.;
    fXoff_edisp = 0.;
    fYoff_edisp = 0.;

    // cut efficiency counter
    fNStats_All = 0;
    fNStats_Rec = 0;
    fNStats_NImagesCut = 0;
    fNStats_Chi2Cut = 0;
    fNStats_CoreErrorCut = 0;
    fNStats_WobbleCut = 0;
    fNStats_WobbleMinCut = 0;
    fNStats_WobbleMaxCut = 0;
}


/*!
  apply cuts on successful reconstruction to input data
*/
bool VTableLookupDataHandler::cut( bool bWrite )
{
    // require at least two images
    if( fNImages < ( int )fTLRunParameter->fTableFillingCut_NImages_min )
    {
        fNStats_NImagesCut++;
        return false;
    }

    // number of reconstructed events
    fNStats_Rec++;

    if( getWobbleOffset() < 0. || getWobbleOffset() > fTLRunParameter->fTableFillingCut_WobbleCut_max )
    {
        // do not cut away single telescope images and if redo steroreconstruction is requested
        if( fNImages != 1 && !fTLRunParameter->fRerunStereoReconstruction )
        {
            fNStats_WobbleCut++;
            return false;
        }
    }

    if( bWrite )
    {
        if( fMCxoff * fMCxoff + fMCyoff * fMCyoff < fMC_distance_to_cameracenter_min * fMC_distance_to_cameracenter_min )
        {
            fNStats_WobbleMinCut++;
            return false;
        }
        if( fMCxoff * fMCxoff + fMCyoff * fMCyoff > fMC_distance_to_cameracenter_max * fMC_distance_to_cameracenter_max )
        {
            fNStats_WobbleMaxCut++;
            return false;
        }
    }

    return true;
}


bool VTableLookupDataHandler::randomSelected()
{
    // random event selection
    if( fSelectRandom > 0. )
    {
        if( fRandom->Uniform() > fSelectRandom )
        {
            return false;
        }
    }
    return true;
}


void VTableLookupDataHandler::setSelectRandom( double iF, int iS )
{
    if( iF > 1. )
    {
        cout << "VTableLookupDataHandler::setSelectRandom error: random selector outside interval [0,1]: " << iF << endl;
        exit( EXIT_FAILURE );
    }

    fSelectRandom = iF;
    fSelectRandomSeed = iS;
    fRandom->SetSeed( iS );
}

/*

   get noise level per telescope type

*/
map< ULong64_t, double > VTableLookupDataHandler::getNoiseLevel_per_TelescopeType()
{
    map< ULong64_t, double > iNSB_telType;
    map< ULong64_t, double > iNSB_telTypeN;
    for( fList_of_Tel_type_iterator = fList_of_Tel_type.begin(); fList_of_Tel_type_iterator != fList_of_Tel_type.end();
            ++fList_of_Tel_type_iterator )
    {
        iNSB_telType[fList_of_Tel_type_iterator->first] = 0.;
        iNSB_telTypeN[fList_of_Tel_type_iterator->first] = 0.;

        for( unsigned int i = 0; i < fNoiseLevel.size(); i++ )
        {
            if( fNoiseLevel[i] > 0. )
            {
                iNSB_telType[fTel_type[i]]  += fNoiseLevel[i];
                iNSB_telTypeN[fTel_type[i]] += 1.;
            }
        }
    }

    // mean value
    map<ULong64_t, double >::iterator iList_of_Tel_type_iterator;
    for( iList_of_Tel_type_iterator = iNSB_telType.begin(); iList_of_Tel_type_iterator != iNSB_telType.end();
            ++iList_of_Tel_type_iterator )
    {

        if( iNSB_telTypeN.find( iList_of_Tel_type_iterator->first ) != iNSB_telTypeN.end()
                && iNSB_telTypeN[iList_of_Tel_type_iterator->first] > 0. )
        {
            iList_of_Tel_type_iterator->second /= iNSB_telTypeN[iList_of_Tel_type_iterator->first];
        }
    }
    return iNSB_telType;
}


/*

   calculates mean noise level over all telescopes with a valid image

   can use current noise level from time dependent pedestal variations

*/
double VTableLookupDataHandler::calculateMeanNoiseLevel( bool bCurrentNoiseLevel )
{
    double z = 0.;
    double m = 0.;

    // time dependent pedestal variations
    if( bCurrentNoiseLevel )
    {
        for( unsigned int i = 0; i < fCurrentNoiseLevel.size(); i++ )
        {
            // is fImgSel_list ntel long or nimages?
            if( fCurrentNoiseLevel[i] > 0. && fImgSel_list[i] )
            {
                m += fCurrentNoiseLevel[i];
                z++;
            }
        }
    }
    // pedestal variations are constant over whole run
    else
    {
        for( unsigned int i = 0; i < fNoiseLevel.size(); i++ )
        {
            if( fNoiseLevel[i] > 0. )
            {
                m += fNoiseLevel[i];
                z++;
            }
        }
    }
    if( z > 0. )
    {
        return m / z;
    }

    return 0.;
}



/*
   get array pointing

   if array pointing does not exist:
   return most probable telescope elevation (majority vote)
*/
double VTableLookupDataHandler::getTelElevation()
{
    // return array pointing elevation
    if( fArrayPointing_Elevation > 1. )
    {
        return fArrayPointing_Elevation;
    }

    // alternative: return most probable telescope elevation
    vector< unsigned int > i_votes( getNTel(), 0 );

    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        for( unsigned int j = 0; j < getNTel(); j++ )
        {
            if( i != j )
            {
                if( TMath::Abs( fTelElevation[i] - fTelElevation[j] ) < 0.2 )
                {
                    i_votes[i]++;
                }
            }
        }
    }
    // get telescope with maximum votes
    unsigned int i_max = 0;
    for( unsigned int i = 0; i < i_votes.size(); i++ )
    {
        if( i_votes[i] > i_votes[i_max] )
        {
            i_max = i;
        }
    }

    cout << "\treading telescope elevation from telescope " << i_max + 1 << " (from vote casting): ";
    cout << fTelElevation[i_max] << " deg" << endl;

    return fTelElevation[i_max];
}

/*
 * get time gradient data vector
 *
 * used for table filling only
*/
double* VTableLookupDataHandler::getTimeGradient( ULong64_t iTelType )
{
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                ftgrad_x_telType[z] = ftgrad_x[i];
                z++;
            }
        }
    }
    return ftgrad_x_telType;
}


/*
 * used for table filling only
 *
 */
double* VTableLookupDataHandler::getDistanceToCore( ULong64_t iTelType, bool iMC )
{
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        fR_telType[i] = 0.;
        fR_telType_MC[i] = 0.;
    }
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                fR_telType[z] = fR[i];
                fR_telType_MC[z] = fR_MC[i];
                z++;
            }
        }
    }
    if( iMC )
    {
        return fR_telType_MC;
    }
    return fR_telType;
}
/*
 * used for table filling only
 *
 */
double* VTableLookupDataHandler::getDistance( ULong64_t iTelType )
{
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        fDistance_telType[i] = 0.;
    }
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                fDistance_telType[z] = fdist[i];
                z++;
            }
        }
    }
    return fDistance_telType;
}

/*
 * used for table filling only
 *
 */
double* VTableLookupDataHandler::getLoss( ULong64_t iTelType )
{
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        fLoss_telType[i] = 0.;
    }
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                fLoss_telType[z] = floss[i];
                z++;
            }
        }
    }
    return fLoss_telType;
}

/*
 * check if this image should be used for the
 * reconstruction steps
 *
 * usually follows the selection from evndisp
 *
 */
bool VTableLookupDataHandler::doImageQualitySelection( unsigned int iTelID )
{
    if( iTelID >= getNTel() )
    {
        return false;
    }
    if( !fImgSel_list[iTelID] )
    {
        return false;
    }
    if( fTLRunParameter->fQualityCutLevel == 0 )
    {
        return true;
    }

    // additional and stricter quality cuts
    if( fTLRunParameter->fQualityCutLevel == 1 )
    {
        // require that centroids are inside the FOV
        if( fdist[iTelID] > fTelFOV[iTelID] * 0.45 )
        {
            return false;
        }
        // set loss for telescope types to 10%
        if( floss[iTelID] > 0.2 )
        {
            return false;
        }
    }

    return true;
}

/*
 * get an array with image size
 *
 * called while reading lookup tables
 *
 * iSelectedImagesOnly = true: use eventdisplay selection
 *
 */
double* VTableLookupDataHandler::getSize( double iSizeCorrection, bool iSelectedImagesOnly, bool iSize2 )
{
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( iSelectedImagesOnly && !doImageQualitySelection( i ) )
        {
            fsizeCorr[i] = -99.;
            continue;
        }
        if( !iSize2 )
        {
            fsizeCorr[i] = fsize[i] * iSizeCorrection;
        }
        else
        {
            fsizeCorr[i] = fsize2[i] * iSizeCorrection;
        }
    }
    return fsizeCorr;
}

/*
 * get an array with image size
 *
 * called while filling lookup tables
 *
 * iSelectedImagesOnly = true: use eventdisplay selection
 *
 */
double* VTableLookupDataHandler::getSize( double iSizeCorrection,  ULong64_t iTelType, bool iSelectedImagesOnly, bool iSize2 )
{
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            if( iSelectedImagesOnly && !doImageQualitySelection( i ) )
            {
                // ignore all telescopes which are not to be analyzed
                if( fTLRunParameter->fTelToAnalyzeData[i]
                        && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
                {
                    fsize_telType[z] = -99.;
                    z++;
                }
                continue;
            }
            if( !iSize2 )
            {
                fsize_telType[z] = fsize[i] * iSizeCorrection;
            }
            else
            {
                fsize_telType[z] = fsize2[i] * iSizeCorrection;
            }
            z++;
        }
    }
    return fsize_telType;
}


/*
 * used for table filling only
 *
 */
double* VTableLookupDataHandler::getWidth( ULong64_t iTelType )
{
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                fwidth_telType[z] = fwidth[i];
                z++;
            }
        }
    }
    return fwidth_telType;
}

/*
 * used for table filling only
 *
 */
double* VTableLookupDataHandler::getLength( ULong64_t iTelType )
{
    unsigned int z = 0;
    for( unsigned int i = 0; i < getNTel(); i++ )
    {
        if( fTel_type[i] == iTelType )
        {
            // ignore all telescopes which are not to be analyzed
            if( fTLRunParameter->fTelToAnalyzeData[i]
                    && fTLRunParameter->fTelToAnalyzeData[i]->fTelToAnalyze )
            {
                flength_telType[z] = flength[i];
                z++;
            }
        }
    }
    return flength_telType;
}

unsigned int VTableLookupDataHandler::getTelType_arraycounter( unsigned int iTelID )
{
    if( iTelID < fTel_type_counter.size() )
    {
        return fTel_type_counter[iTelID];
    }

    return 999999;
}

/*
 * initialize vector which assigns for each telescope the telescope type counter
 *
 */
void VTableLookupDataHandler::initializeTelTypeVector()
{
    fTel_type_counter.assign( fNTel, 9999 );
    for( unsigned int iTelID = 0; iTelID < fNTel; iTelID++ )
    {
        unsigned int z = 0;
        for( fList_of_Tel_type_iterator = fList_of_Tel_type.begin();
                fList_of_Tel_type_iterator != fList_of_Tel_type.end(); fList_of_Tel_type_iterator++ )
        {
            if( fTel_type[iTelID] == fList_of_Tel_type_iterator->first )
            {
                fTel_type_counter[iTelID] = z;
            }
            z++;
        }
    }
}
