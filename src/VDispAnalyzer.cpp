/* \class VDispAnalyzer


     calculate displacement as function of width, length, size, ...

     call different methods to get disp (lookup tables, neural network, boosted decision trees)

     allows also to reconstruct the energy, core, dispError (implemented for TMVA disp only)

*/

#include "VDispAnalyzer.h"

VDispAnalyzer::VDispAnalyzer()
{
    fDispTableAnalyzer = 0;
    fTMVADispAnalyzer = 0;

    fDispMethod = "NOMETHODDEFINED";

    f_disp = -999.;
    f_dispE = -999.;
    f_dispDiff = -999.;
    f_xs = -999.;
    f_ys = -999.;
    f_angdiff = 0.;

    fdisp_energy = -9999.;
    fdisp_energy_chi = -9999.;
    fdisp_energy_dEs = -9999.;
    fdisp_energy_median = -9999.;
    fdisp_energy_medianAbsoluteError = -9999.;
    fdisp_energy_NT = 0;
    fdisp_energyQL = -1;
    fdisp_sum_abs_weigth = 0.;

    setQualityCuts();
    setDispErrorWeighting();
    setDebug( false );
}
void VDispAnalyzer::setTelescopeTypeList( vector<ULong64_t> iTelescopeTypeList )
{
    fTelescopeTypeList = iTelescopeTypeList;
}

/*
 * initialize the disp methods
 *
 * e.g. read lookup tables or weight files for BDTs or NN
 *
 */
bool VDispAnalyzer::initialize( string iFile, string iDispMethod, string iDispType )
{
    fDispMethod = iDispMethod;

    // disp is calculated using simple lookup tables
    if( fDispMethod == "DISPTABLES" )
    {
        fDispTableAnalyzer = new VDispTableAnalyzer( iFile );
        if( fDispTableAnalyzer->isZombie() )
        {
            setZombie( true );
        }
        else
        {
            setZombie( false );
        }
    }
    // disp is calculated using TMVA BDTs
    else if( fDispMethod == "TMVABDT" )
    {
        fTMVADispAnalyzer = new VTMVADispAnalyzer( iFile, fTelescopeTypeList, iDispType );
        if( fTMVADispAnalyzer->isZombie() )
        {
            setZombie( true );
        }
        else
        {
            setZombie( false );
        }
    }
    else
    {
        cout << "VDispAnalyzer::initialize ERROR: unknown disp method: " << fDispMethod << endl;
        return false;
    }

    if( isZombie() )
    {
        cout << "VDispAnalyzer::initialize ERROR initializing method " << fDispMethod << endl;
        cout << "\t no input file: " << iFile << endl;
        cout << "\t exiting..." << endl;
        exit( EXIT_FAILURE );
    }

    return true;
}

/*
 * finish orderly
 *
 *
 */
void VDispAnalyzer::terminate()
{
    if( fDispTableAnalyzer )
    {
        fDispTableAnalyzer->terminate();
    }

    if( fTMVADispAnalyzer )
    {
        fTMVADispAnalyzer->terminate();
    }
}

/*
 * calculate disp from the parameters given
 *
 * (for one individual image)
 *
*/
float VDispAnalyzer::evaluate( float iWidth, float iLength, float iAsymm, float iDist, float iSize,
                               float iPedvar, float iTGrad, float iLoss, float icen_x, float icen_y,
                               float xoff_4, float yoff_4, ULong64_t iTelType,
                               float iZe, float iAz, float iRcore, float iFui, float iNtubes, bool b2D )
{
    f_disp = -99.;

    if( fDispTableAnalyzer )
    {
        f_disp = fDispTableAnalyzer->evaluate( iWidth, iLength, iSize, iPedvar, iZe, iAz, b2D );
    }
    else if( fTMVADispAnalyzer )
    {
        f_disp = fTMVADispAnalyzer->evaluate( iWidth, iLength, iSize, iAsymm, iLoss, iTGrad,
                                              icen_x, icen_y, xoff_4, yoff_4, iTelType,
                                              iZe, iAz, iRcore, -1., iDist, iFui, iNtubes );
    }
    return f_disp;

}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/*
 * return permutation vector with all possible pairs of signs
 *
 */
vector< vector< float > > VDispAnalyzer::get_sign_permutation_vector( unsigned int x_size )
{
    unsigned int ncombinations = ( 1 << x_size ) - 1;
    vector< vector< float > > i_sign;
    for( unsigned int i = 0; i <= ncombinations; i++ )
    {
        int y = i;
        int i_counter = 0;
        vector< float > i_temp_v( x_size, 1. );
        while( y > 0 )
        {
            if( ( y & 1 ) == 1 )
            {
                i_temp_v[i_counter] = -1.;
            }
            y = y >> 1;
            i_counter++;
        }
        i_sign.push_back( i_temp_v );
    }
    return i_sign;
}

/*
 * calculate smallest diff element of a set of disp
 * values
 *
 */
unsigned int VDispAnalyzer::find_smallest_diff_element(
    vector< vector< float > > i_sign,
    vector< float > x, vector< float > y,
    vector< float > cosphi, vector< float > sinphi,
    vector< float > v_disp, vector< float > v_weight )
{
    vector< float > v_disp_diff( i_sign.size(), 0. );
    vector< float > v_dist( i_sign.size(), 0. );
    for( unsigned int s = 0; s < i_sign.size(); s++ )
    {
        vector< float > v_xs;
        vector< float > v_ys;
        float xs = 0.;
        float ys = 0.;
        float disp_diff = 0.;
        for( unsigned int i = 0; i < x.size(); i++ )
        {
            v_xs.push_back( x[i] - i_sign[s][i] * v_disp[i] * cosphi[i] );
            v_ys.push_back( y[i] - i_sign[s][i] * v_disp[i] * sinphi[i] );
        }
        calculateMeanShowerDirection( v_xs, v_ys, v_weight, xs, ys, disp_diff, v_xs.size() );
        v_disp_diff[s] = disp_diff;
        v_dist[s] = sqrt( xs * xs + ys * ys );
    }
    unsigned int i_mean_element = 0;
    float i_smallest_dist = 1.e20;

    for( unsigned int s = 0; s < v_disp_diff.size(); s++ )
    {
        if( v_disp_diff[s] < i_smallest_dist )
        {
            i_mean_element = s;
            i_smallest_dist = v_disp_diff[s];
        }
    }
    return i_mean_element;
}

/*
 * return index of N largest values of a vector
 *
 */
vector< unsigned int > VDispAnalyzer::get_largest_weight_index( vector< float > w, unsigned int N )
{
    multimap< float, unsigned int, greater<float> > w_sorted;
    for( unsigned int i = 0; i < w.size(); i++ )
    {
        w_sorted.insert( pair< float, unsigned int>( w[i], i ) );
    }
    vector< unsigned int > rev_sorted_index;
    for( multimap<float, unsigned int>::iterator it = w_sorted.begin(); it != w_sorted.end(); ++it )
    {
        rev_sorted_index.push_back( it->second );
        if( rev_sorted_index.size() > N )
        {
            break;
        }
    }
    return rev_sorted_index;
}


/*
 * calculate direction coordinates (x,y) from an array of points using
 * disp and centroids
 *
 * return values:
 * xs, ys:           reconstructed shower direction
 * dispdiff:
 *
 * input values:
 * x, y:             centroids of images
 * cosphi, shinphi:  orientation of image axes
 * v_disp:           displacement [deg]
 * v_weigth:         weight using in calculation of average direction
 *
 *
 */
void VDispAnalyzer::calculateMeanDirection( float& xs, float& ys,
        vector< float > x, vector< float > y,
        vector< float > cosphi, vector< float > sinphi,
        vector< float > v_disp, vector< float > v_weight,
        float& dispdiff,
        float x_off4, float y_off4 )
{
    xs = -99999.;
    ys = -99999.;
    dispdiff = -9999.;

    ////////////////////////////////////////////////////
    // note: first part (calculating smallest difference between solutions)
    // works only for a maximum of 4 images
    // for more than 4 images: seed value from different reconstruction method
    // should be given.
    if( x.size() > 4 && x_off4 < -998. && y_off4 < -998. )
    {
        return;
    }

    //////////////////////////////////////////////////////////
    // calculate (average) angle between the image lines for
    f_angdiff = 0.;
    float fmean_iangdiffN = 0.;

    if( cosphi.size() > 1 )
    {
        // calculate average angle between image lines
        for( unsigned int ii = 0; ii < sinphi.size(); ii++ )
        {
            for( unsigned int jj = 1; jj < sinphi.size(); jj++ )
            {
                if( ii >= jj )
                {
                    continue;
                }
                float mi = sinphi[ii] / cosphi[ii];
                float mj = sinphi[jj] / cosphi[jj];
                if( isinf( mi ) || isinf( mj ) )
                {
                    continue;
                }
                float iangdiff = fabs( atan( mj ) - atan( mi ) );
                // mean angle between images
                if( iangdiff < 90. * TMath::DegToRad() )
                {
                    f_angdiff += iangdiff * TMath::RadToDeg();
                }
                else
                {
                    f_angdiff += ( 180. - iangdiff * TMath::RadToDeg() );
                }
                fmean_iangdiffN++;
            }
        }
        if( fmean_iangdiffN > 0. )
        {
            f_angdiff /= fmean_iangdiffN;
        }
    }
    // multiplicity 1 (fix to artificially large value)
    else
    {
        fmean_iangdiffN = 1.;
        f_angdiff = 180.;
    }
    // check for close to parallel lines for image multiplicity 2
    // (not so important for disp direction,
    //  but note that core reconstruction
    //  is still done the conventional way)
    // Using LIN51 leads to an increase in effective area at large
    // energies and large off-axis angle, but significantly worse
    // angular resolution
    // LIN51 if( f_angdiff < fAxesAngles_min && fmean_iangdiffN < 2.01 )
    if( f_angdiff < fAxesAngles_min )
    {
        return;
    }

    ///////////////////////////////////////////////////////////////////////
    // method: disp table
    ///////////////////////////////////////////////////////////////////////
    if( fDispTableAnalyzer )
    {
        fDispTableAnalyzer->calculateMeanDirection( xs, ys, x, y, cosphi, sinphi, v_disp, v_weight );
    }
    ///////////////////////////////////////////////////////////////////////
    // method: BDTs
    ///////////////////////////////////////////////////////////////////////
    // head/tail uncertainty
    // search for combination of images with smallest differences
    // in reconstructed images
    fdisp_xs_T.clear();
    fdisp_ys_T.clear();
    fdisp_xs_T.assign( v_weight.size(), 0. );
    fdisp_ys_T.assign( v_weight.size(), 0. );

    unsigned int i_smallest_diff_element = 0;
    // search for cloud with smallest dist
    // 5 images and less: use all disp solutions
    // more images: use 5 images with largest weight
    if( x.size() < 6 )
    {
        vector< vector< float > > i_sign = get_sign_permutation_vector( x.size() );
        i_smallest_diff_element = find_smallest_diff_element(
                                      i_sign, x, y, cosphi, sinphi, v_disp, v_weight );

        for( unsigned int ii = 0; ii < x.size(); ii++ )
        {
            fdisp_xs_T[ii] = x[ii] - i_sign[i_smallest_diff_element][ii] * v_disp[ii] * cosphi[ii];
            fdisp_ys_T[ii] = y[ii] - i_sign[i_smallest_diff_element][ii] * v_disp[ii] * sinphi[ii];
        }
        calculateMeanShowerDirection( fdisp_xs_T, fdisp_ys_T, v_weight, xs, ys, dispdiff, fdisp_xs_T.size() );
    }
    // more images
    // (MVA disp analyzer using
    // stereo reconstruction as starting
    // value
    else if( x_off4 > -998. && y_off4 > -998. )
    {
        // vector with disp directions per image
        fdisp_xs_T.clear();
        fdisp_ys_T.clear();
        fdisp_xs_T.assign( v_weight.size(), 0. );
        fdisp_ys_T.assign( v_weight.size(), 0. );

        float x1 = 0.;
        float x2 = 0.;
        float y1 = 0.;
        float y2 = 0.;

        for( unsigned int ii = 0; ii < v_weight.size(); ii++ )
        {
            x1 = x[ii] - v_disp[ii] * cosphi[ii];
            x2 = x[ii] + v_disp[ii] * cosphi[ii];

            y1 = y[ii] - v_disp[ii] * sinphi[ii];
            y2 = y[ii] + v_disp[ii] * sinphi[ii];

            // check solution closest to starting value
            // (should work properly here, as these are
            // all events with multiplicity > 4)
            if( sqrt( ( x1 - x_off4 ) * ( x1 - x_off4 ) + ( y1 + y_off4 ) * ( y1 + y_off4 ) )
                    < sqrt( ( x2 - x_off4 ) * ( x2 - x_off4 ) + ( y2 + y_off4 ) * ( y2 + y_off4 ) ) )
            {
                fdisp_xs_T[ii] = x1;
                fdisp_ys_T[ii] = y1;
            }
            else
            {
                fdisp_xs_T[ii] = x2;
                fdisp_ys_T[ii] = y2;
                v_disp[ii] = -1.*TMath::Abs( v_disp[ii] );
            }
        }
        calculateMeanShowerDirection( fdisp_xs_T, fdisp_ys_T, v_weight, xs, ys, dispdiff, fdisp_xs_T.size() );
    }

    // apply a completely unnecessary sign flip
    if( ys > -9998. )
    {
        ys *= -1.;
    }
}

/*
 * calculate average shower directions from a vector of disp direction
 *
 * use given weights
 *
 * (internal function)
 *
*/
void VDispAnalyzer::calculateMeanShowerDirection( vector< float > v_x, vector< float > v_y, vector< float > v_weight,
        float& xs, float& ys, float& dispdiff,
        unsigned int iMaxN )
{
    xs = 0.;
    ys = 0.;
    dispdiff = 0.;
    float d_w_sum = 0.;

    if( iMaxN > v_x.size() )
    {
        cout << "VDispAnalyzer::calculateMeanShowerDirection error: ";
        cout << "invalid vector size " << endl;
        exit( EXIT_FAILURE );
    }

    // single image
    if( iMaxN == 1 && v_x.size() == 1 && v_y.size() == 1 && v_weight.size() == 1 )
    {
        dispdiff = 0.;
        xs = v_x[0];
        ys = v_y[0];
        return;
    }

    float z = 0.;
    for( unsigned int n = 0; n < iMaxN; n++ )
    {
        for( unsigned int m = n + 1; m < iMaxN; m++ )
        {
            dispdiff += sqrt( ( v_x[n] - v_x[m] ) * ( v_x[n] - v_x[m] )
                              + ( v_y[n] - v_y[m] ) * ( v_y[n] - v_y[m] ) )
                        * v_weight[n] * v_weight[m];
            z += v_weight[n] * v_weight[m];
        }
        xs += v_x[n] * v_weight[n];
        ys += v_y[n] * v_weight[n];
        d_w_sum += v_weight[n];
    }
    if( d_w_sum > 0. )
    {
        xs /= d_w_sum;
        ys /= d_w_sum;
    }
    if( z > 0. )
    {
        dispdiff /= z;
    }
    if( v_x.size() < 2 && iMaxN >= 2 )
    {
        dispdiff = -9999.;
        xs = -99999.;
        ys = -99999.;
    }
}

/*
 * calculate mean disp direction using as input C arrays
 * (used primarily from lookup table code)
 *
 * called from VTableLookupDataHandler::doStereoReconstruction()
 *
 */
void VDispAnalyzer::calculateMeanDirection( unsigned int i_ntel,
        float iArrayElevation,
        float iArrayAzimuth,
        ULong64_t* iTelType,
        double* img_size,
        double* img_cen_x,
        double* img_cen_y,
        double* img_cosphi,
        double* img_sinphi,
        double* img_width,
        double* img_length,
        double* img_asym,
        double* img_tgrad,
        double* img_loss,
        int* img_ntubes,
        double* img_weight,
        double xoff_4,
        double yoff_4,
        vector< float > dispErrorT,
        double* img_fui )
{
    // reset values from previous event
    f_disp = -99.;
    f_dispDiff = -99.;
    f_xs = -99.;
    f_ys = -99.;

    // make sure that all data arrays exist
    if( !img_size || !img_cen_x || !img_cen_y
            || !img_cosphi || !img_sinphi
            || !img_width || !img_length
            || !img_asym || !img_tgrad
            || !img_loss || !img_ntubes
            || !img_weight || !img_fui )
    {
        return;
    }

    float disp = 0.;
    vector< float > v_disp;
    vector< unsigned int > v_displist;
    vector< float > v_weight;
    vector< float > x;
    vector< float > y;
    vector< float > cosphi;
    vector< float > sinphi;

    //////////////////////////////
    // loop over all telescopes and calculate disp per telescope
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        // quality cuts
        if( img_size[i] > 0. && img_length[i] > 0.
                && sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistance_max
                && img_loss[i] < floss_max
                && img_fui[i] > fFui_min
                && ( fdistanceQC_max == 0 || sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistanceQC_max[i] ) )
        {
            disp = evaluate( ( float )img_width[i], ( float )img_length[i], ( float )img_asym[i],
                             ( float )sqrt( img_cen_x[i] * img_cen_x[i] + img_cen_y[i] * img_cen_y[i] ),
                             ( float )img_size[i], 1., ( float )img_tgrad[i], ( float )img_loss[i],
                             ( float )img_cen_x[i], ( float )img_cen_y[i],
                             ( float )xoff_4, ( float )yoff_4, iTelType[i],
                             ( float )( 90. - iArrayElevation ), ( float )iArrayAzimuth,
                             -99., ( float )img_fui[i], ( float )img_ntubes[i] );
            if( disp < -98. )
            {
                continue;
            }
            // use time gradient to get right directory for image
            if( img_tgrad[i] < 0. )
            {
                disp *= -1.;
            }
            v_disp.push_back( disp );
            v_displist.push_back( i );

            // use estimated uncertainty on disp direction reconstruction as
            // weight: exponential coefficient ad-hoc, not a result of
            // optimisation
            if( fDispErrorWeighting )
            {
                if( i < dispErrorT.size() )
                {
                    v_weight.push_back( exp( -1. * fDispErrorExponential * TMath::Abs( dispErrorT[i] ) ) );
                }
                else
                {
                    v_weight.push_back( -999. );
                }
            }
            // no estimated uncertainties: use image size and geometry
            else
            {
                v_weight.push_back( img_size[i] * img_weight[i] * ( 1. - img_width[i] / img_length[i] )
                                    * img_size[i] * img_weight[i] * ( 1. - img_width[i] / img_length[i] )
                                    * img_fui[i] * img_fui[i] * img_fui[i] * img_fui[i] );
            }

            x.push_back( img_cen_x[i] );
            y.push_back( img_cen_y[i] );
            cosphi.push_back( img_cosphi[i] );
            sinphi.push_back( img_sinphi[i] );
        }
    }

    // calculate expected direction
    calculateMeanDirection( f_xs, f_ys,
                            x, y, cosphi, sinphi,
                            v_disp, v_weight,
                            f_dispDiff, xoff_4, yoff_4 );
    fdisp_xy_weight_T = v_weight;
    fdisp_sum_abs_weigth = 0.;
    for( unsigned int i = 0; i < fdisp_xy_weight_T.size(); i++ )
    {
        fdisp_sum_abs_weigth += TMath::Abs( fdisp_xy_weight_T[i] );
    }
    fdisp_T = v_disp;
    fdisplist_T = v_displist;
}

float VDispAnalyzer::getDispErrorT( unsigned int iTelescopeNumber )
{
    if( iTelescopeNumber < fdisp_error_T.size() && fdisp_error_T[iTelescopeNumber] > 0. )
    {
        return fdisp_error_T[iTelescopeNumber];
    }

    return -9999.;
}


/*
 * calculate expected error for angular reconstruction
 * using the disp method
 *
 * one value per telescope
 *
 * called from VTableLookupDataHandler::doStereoReconstruction()
 *
 */
void VDispAnalyzer::calculateExpectedDirectionError( unsigned int i_ntel,
        float iArrayElevation,
        float iArrayAzimuth,
        ULong64_t* iTelType,
        double* img_size,
        double* img_cen_x,
        double* img_cen_y,
        double* img_cosphi,
        double* img_sinphi,
        double* img_width,
        double* img_length,
        double* img_asym,
        double* img_tgrad,
        double* img_loss,
        int* img_ntubes,
        double* img_weight,
        double xoff_4,
        double yoff_4,
        double* img_fui )
{
    // make sure that all data arrays exist
    if( !img_size || !img_cen_x || !img_cen_y
            || !img_cosphi || !img_sinphi
            || !img_width || !img_length
            || !img_asym || !img_tgrad
            || !img_loss || !img_ntubes
            || !img_weight || !img_fui )
    {
        return;
    }
    fdisp_error_T.clear();
    fdisp_error_T.assign( i_ntel, -99. );

    //////////////////////////////
    // loop over all telescopes and calculate disp per telescope
    for( unsigned int i = 0; i < i_ntel; i++ )
    {

        // quality cuts
        if( img_size[i] > 0. && img_length[i] > 0.
                && sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistance_max
                && img_loss[i] < floss_max
                && img_fui[i] > fFui_min
                && ( fdistanceQC_max == 0 || sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistanceQC_max[i] ) )
        {
            fdisp_error_T[i] = evaluate( ( float )img_width[i], ( float )img_length[i], ( float )img_asym[i],
                                         ( float )sqrt( img_cen_x[i] * img_cen_x[i] + img_cen_y[i] * img_cen_y[i] ),
                                         ( float )img_size[i], 1., ( float )img_tgrad[i], ( float )img_loss[i],
                                         ( float )img_cen_x[i], ( float )img_cen_y[i],
                                         ( float )xoff_4, ( float )yoff_4, iTelType[i],
                                         ( float )( 90. - iArrayElevation ), ( float )iArrayAzimuth,
                                         -99., ( float )img_fui[i], ( float )img_ntubes[i] );
        }
    }
}


/*
 * calculate x coordinate from disp, centroid, and image line orientation
 *
*/
float VDispAnalyzer::getXcoordinate_disp( unsigned int ii )
{
    // disp table analysis
    if( fDispTableAnalyzer )
    {
        return fDispTableAnalyzer->getXcoordinate_disp( ii );
    }
    // other disp analyses
    else if( ii < fdisp_xs_T.size() )
    {
        return fdisp_xs_T[ii];
    }

    return -99.;
}

/*
 * calculate y coordinate from disp, centroid, and image line orientation
 *
*/
float VDispAnalyzer::getYcoordinate_disp( unsigned int ii )
{
    if( fDispTableAnalyzer )
    {
        return fDispTableAnalyzer->getYcoordinate_disp( ii );
    }
    // other disp analyses
    // (include sign flip convention)
    else if( ii < fdisp_ys_T.size() )
    {
        return -1.*fdisp_ys_T[ii];
    }

    return -99.;
}

/**************************************************************************************
 * everything below is related to dispENERGY
 **************************************************************************************/

/*
 * calculate average energy for the given event
 * using the dispEnergy algorithms
 *
 */
void VDispAnalyzer::calculateEnergies( unsigned int i_ntel,
                                       float iArrayElevation,
                                       float iArrayAzimuth,
                                       ULong64_t* iTelType,
                                       double* img_size,
                                       double* img_cen_x,
                                       double* img_cen_y,
                                       double* img_cosphi,
                                       double* img_sinphi,
                                       double* img_width,
                                       double* img_length,
                                       double* img_asym,
                                       double* img_tgrad,
                                       double* img_loss,
                                       int* img_ntubes,
                                       double* img_weight,
                                       double xoff_4,
                                       double yoff_4,
                                       double* iRcore,
                                       double iEHeight,
                                       double iMCEnergy,
                                       double* img_fui )
{
    fdisp_energy = -99.;
    fdisp_energy_chi = 1.e-10;
    fdisp_energy_dEs = 1.e-10;
    fdisp_energy_T.clear();
    fdisp_energy_T.assign( i_ntel, -99. );
    fdisp_energy_NT = 0;
    fdisp_energyQL = -1;

    // make sure that all data arrays exist
    if( !img_size || !img_cen_x || !img_cen_y
            || !img_cosphi || !img_sinphi
            || !img_width || !img_length
            || !img_asym || !img_tgrad
            || !img_loss || !img_ntubes
            || !img_weight || !iRcore
            || !img_fui )
    {
        return;
    }

    ////////////////////////////////////////////
    // calculate for each image an energy

    // counter for good energy values
    float z = 0.;
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( img_size[i] > 0. && iRcore[i] > 0. && iArrayElevation > 0.
                && sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistance_max
                && img_loss[i] < floss_max
                && img_fui[i] > fFui_min
                && ( fdistanceQC_max == 0 || sqrt( img_cen_x[i]*img_cen_x[i] + img_cen_y[i]*img_cen_y[i] ) < fdistanceQC_max[i] ) )
        {
            fdisp_energy_T[i] = fTMVADispAnalyzer->evaluate(
                                    ( float )img_width[i], ( float )img_length[i],
                                    ( float )img_size[i], ( float )img_asym[i],
                                    ( float )img_loss[i], ( float )img_tgrad[i],
                                    ( float )img_cen_x[i], ( float )img_cen_y[i],
                                    ( float )xoff_4, ( float )yoff_4, ( ULong64_t )iTelType[i],
                                    ( float )( 90. - iArrayElevation ), ( float )iArrayAzimuth,
                                    ( float )iRcore[i], ( float )iEHeight,
                                    ( float )sqrt( img_cen_x[i] * img_cen_x[i] + img_cen_y[i] * img_cen_y[i] ),
                                    ( float )img_fui[i], ( float )img_ntubes[i] );

            if( fDebug )
            {
                cout << "VDispAnalyzer::calculateEnergies: tel " << i << " (teltype " << ( ULong64_t )iTelType[i] << ") ";
                cout << "  size " << img_size[i] << " R " << iRcore[i] << "\t loss " << img_loss[i] << endl;
                cout << "\t\t energy : " << fdisp_energy_T[i] << "\t Erec/MCe0: " << fdisp_energy_T[i] / iMCEnergy << endl;
            }
            z++;
        }
        else
        {
            fdisp_energy_T[i] = -99.;
        }
    }
    // check that there is a minimum number
    // of good energy values
    if( z < 1.e-4 )
    {
        return;
    }

    ///////////////////
    // calculate average energy
    // (require at least two telescopes)

    // fill a 'clean' vector with good energy per telescopes
    vector< double > energy_tel;
    vector< double > energy_weight;
    vector< double > iR;
    vector< double > iS;
    vector< double > iW;
    vector< double > iL;
    vector< double > iT;
    vector< double > iLe;
    vector< double > iWi;

    for( unsigned int i = 0; i < fdisp_energy_T.size(); i++ )
    {
        if( fdisp_energy_T[i] > 0. && img_weight[i] > 0. && img_size[i] > 0. )
        {
            energy_tel.push_back( fdisp_energy_T[i] );
            // Weighting with size leads to a long tail towards large eres/mce0
            //              energy_weight.push_back( img_size[i] * img_weight[i] * img_size[i] * img_weight[i] );
            energy_weight.push_back( img_weight[i] * img_weight[i] );
            if( fDebug )
            {
                iR.push_back( iRcore[i] );
                iW.push_back( img_weight[i] );
                iS.push_back( img_size[i] );
                iL.push_back( img_loss[i] );
                iT.push_back( img_tgrad[i] );
                iLe.push_back( img_length[i] );
                iWi.push_back( img_width[i] );
            }
        }
    }

    // Occasionally one energy is significantly off and distorts the mean.
    // therefore: get rid of N sigma outliers
    // use robust statistics (median and median absolute error)
    // Note: applied only to larger events > 4 telescopes
    fdisp_energy_median = TMath::Median( energy_tel.size(), &energy_tel[0] );
    fdisp_energy_medianAbsoluteError = VStatistics::getMedianAbsoluteError( energy_tel, fdisp_energy_median );
    double w = 0.;
    unsigned int n2 = 0;
    fdisp_energy = 0.;
    for( unsigned int j = 0; j < energy_tel.size(); j++ )
    {
        if( energy_tel.size() < 5
                || TMath::Abs( energy_tel[j] - fdisp_energy_median ) < fdisp_energy_medianAbsoluteError * 3. )
        {
            fdisp_energy += energy_tel[j] * energy_weight[j];
            w += energy_weight[j];
            n2++;
        }
    }
    fdisp_energy_NT = energy_tel.size();
    // check minimum number of valid energies
    /*     if( energy_tel.size() <= 4 && w > 0. )
         {
             fdisp_energy /= w;
             fdisp_energy_NT = n2;
         }  */
    // use median for energy estimation (removes outliers)
    // (for all cases)
    if( n2 >= fNImages_min )
    {
        fdisp_energy = fdisp_energy_median;
        if( n2 == 1 )
        {
            fdisp_energyQL = 1;
        }
        else
        {
            fdisp_energyQL = 0;
        }
    }
    else
    {
        fdisp_energy = -99.;
        fdisp_energy_NT = 0;
        fdisp_energyQL = -1;
    }

    ///////////////////
    // calculate chi2 and dE
    // (note: different definition for dE
    //        than in lookup table code)
    z = 0.;
    fdisp_energy_chi = 0.;
    fdisp_energy_dEs = 0.;
    for( unsigned int i = 0; i < fdisp_energy_T.size(); i++ )
    {
        if( fdisp_energy_T[i] > 0. && fdisp_energy > 0. )
        {
            fdisp_energy_chi += ( fdisp_energy_T[i] - fdisp_energy ) *
                                ( fdisp_energy_T[i] - fdisp_energy ) /
                                fdisp_energy / fdisp_energy;
            fdisp_energy_dEs += TMath::Abs( fdisp_energy_T[i] - fdisp_energy ) / fdisp_energy;
            z++;
        }
    }
    if( z > 1.5 )
    {
        fdisp_energy_chi /= ( z - 1. );
        fdisp_energy_dEs /= z;
    }
    else
    {
        fdisp_energy_chi = 1.e-10;
        fdisp_energy_dEs = 1.e-10;
    }
}

///////////////////////////////////////////////////////////////

float VDispAnalyzer::getEnergy()
{
    return fdisp_energy;
}

float VDispAnalyzer::getEnergyChi2()
{
    return fdisp_energy_chi;
}

float VDispAnalyzer::getEnergydES()
{
    return fdisp_energy_dEs;
}

float VDispAnalyzer::getEnergyMedian()
{
    return fdisp_energy_median;
}

float VDispAnalyzer::getEnergyMedianAbsoluteError()
{
    return fdisp_energy_medianAbsoluteError;
}

float VDispAnalyzer::getEnergyT( unsigned int iTelescopeNumber )
{
    if( iTelescopeNumber < fdisp_energy_T.size() )
    {
        return fdisp_energy_T[iTelescopeNumber];
    }

    return -9999.;
}

/**************************************************************************************
 * everything below is related to dispCore
 **************************************************************************************/

void VDispAnalyzer::calculateCore( unsigned int i_ntel,
                                   float iArrayElevation,
                                   float iArrayAzimuth,
                                   double* iTelX,
                                   double* iTelY,
                                   double* iTelZ,
                                   ULong64_t* iTelType,
                                   double* img_size,
                                   double* img_cen_x,
                                   double* img_cen_y,
                                   double* img_cosphi,
                                   double* img_sinphi,
                                   double* img_width,
                                   double* img_length,
                                   double* img_asym,
                                   double* img_tgrad,
                                   double* img_loss,
                                   int* img_ntubes,
                                   double* img_weight,
                                   double xoff_4,
                                   double yoff_4,
                                   double* iRcore,
                                   double xcoreSR,
                                   double ycoreSR,
                                   double xs,
                                   double ys,
                                   double* img_fui )
{
    fdisp_core_T.clear();
    fdisp_core_T.assign( i_ntel, -99. );

    // make sure that all data arrays exist
    if( !img_size || !img_cen_x || !img_cen_y
            || !img_cosphi || !img_sinphi
            || !img_width || !img_length
            || !img_asym || !img_tgrad
            || !img_loss || !img_ntubes
            || !img_weight || !iRcore
            || !iTelX || !iTelY || !iTelZ || !img_fui )
    {
        return;
    }

    ////////////////////////////////////////////
    // calculate for each image the core distance

    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( img_size[i] > 0. && iArrayElevation > 0. )
        {
            fdisp_core_T[i] = fTMVADispAnalyzer->evaluate(
                                  ( float )img_width[i], ( float )img_length[i],
                                  ( float )img_size[i], ( float )img_asym[i],
                                  ( float )img_loss[i], ( float )img_tgrad[i],
                                  ( float )img_cen_x[i], ( float )img_cen_y[i],
                                  ( float )xoff_4, ( float )yoff_4, ( ULong64_t )iTelType[i],
                                  ( float )( 90. - iArrayElevation ), ( float )iArrayAzimuth,
                                  ( float )iRcore[i], -1.,
                                  ( float )sqrt( img_cen_x[i] * img_cen_x[i] + img_cen_y[i] * img_cen_y[i] ),
                                  ( float )img_fui[i], ( float )img_ntubes[i] );
        }
        else
        {
            fdisp_core_T[i] = -99.;
        }
    }

    return;

    ///////////////////////////////////////////////
    // core position in shower coordinates

    VGrIsuAnalyzer i_A;

    // calculated telescope positions in shower coordinates
    float i_xcos = sin( ( 90. - iArrayElevation ) / TMath::RadToDeg() ) * sin( ( iArrayAzimuth - 180. ) / TMath::RadToDeg() );
    float i_ycos = sin( ( 90. - iArrayElevation ) / TMath::RadToDeg() ) * cos( ( iArrayAzimuth - 180. ) / TMath::RadToDeg() );
    float i_xrot, i_yrot, i_zrot = 0.;
    float xcoreSC = 0.;
    float ycoreSC = 0.;
    i_A.tel_impact( i_xcos, i_ycos, xcoreSR, ycoreSR, 0., &xcoreSC, &ycoreSC, &i_zrot, false );

    float m = 0.;
    float theta = 0.;
    float x = 0.;
    float y = 0.;
    float x1 = 0.;
    float x2 = 0.;
    float y1 = 0.;
    float y2 = 0.;

    float iweight = 0.;
    float xw = 0.;
    float yw = 0.;
    float ww = 0.;

    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( fdisp_core_T[i] > 0. && img_cen_x[i] - xs > 0. && xcoreSC > -9000. && ycoreSC > -9000. )
        {
            // telescope coordinates
            // shower coordinates (telecope pointing)
            i_A.tel_impact( i_xcos, i_ycos, iTelX[i], iTelY[i], iTelZ[i], &i_xrot, &i_yrot, &i_zrot, false );

            m = -1. * ( img_cen_y[i] - ys ) / ( img_cen_x[i] - xs );

            theta = atan( m );

            // not clear which sign to use
            // simple approximation: assume
            // that intersection by line solution
            // is not hugely wrong; look for closest
            // solution
            x1 = i_xrot + fdisp_core_T[i] * cos( theta );
            y1 = i_yrot + fdisp_core_T[i] * sin( theta );
            x2 = i_xrot - fdisp_core_T[i] * cos( theta );
            y2 = i_yrot - fdisp_core_T[i] * sin( theta );
            if( sqrt( ( x1 - xcoreSC ) * ( x1 - xcoreSC )
                      + ( y1 - ycoreSC ) * ( y1 - ycoreSC ) )
                    < sqrt( ( x2 - xcoreSC ) * ( x2 - xcoreSC )
                            + ( y2 - ycoreSC ) * ( y2 - ycoreSC ) ) )
            {
                x = x1;
                y = y1;
            }
            else
            {
                x = x2;
                y = y2;
            }

            // weight
            // (not clear if sum is a good choice,
            //  as we do not do a dc-to-pe calibration)
            iweight = img_size[i];
            iweight *= ( 1. - img_width[i] / img_length[i] );
            iweight *= iweight;

            xw += x * iweight;
            yw += y * iweight;
            ww += iweight;
        }
    }

    // average core position (in shower coordinates)
    if( ww > 0. )
    {
        xw /= ww;
        yw /= ww;
    }

    // calculate core distance for each telescopes
    // (simple - as we are working in shower coordinates)
    for( unsigned int i = 0; i < i_ntel; i++ )
    {
        if( fdisp_core_T[i] > 0. && ww > 0. && img_cen_x[i] > 0. )
        {
            i_A.tel_impact( i_xcos, i_ycos, iTelX[i], iTelY[i], iTelZ[i], &i_xrot, &i_yrot, &i_zrot, false );
            fdisp_core_T[i] = sqrt( ( xw - i_xrot ) * ( xw - i_xrot )
                                    + ( yw - i_yrot ) * ( yw - i_yrot ) );
        }
    }

}

float VDispAnalyzer::getCoreDistance( unsigned int iTelescopeNumber )
{
    if( iTelescopeNumber < fdisp_core_T.size() )
    {
        return fdisp_core_T[iTelescopeNumber];
    }
    return -99.;
}
