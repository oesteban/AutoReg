/*
 Copyright (c) 2010, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM)
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of the BIT-UPM, nor the names of its contributors
 may be used to endorse or promote products derived from this software
 without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY Oscar Esteban ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL OSCAR ESTEBAN BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "autoreg.h"

int main(int argc, char *argv[]) {
    typedef MIRegistration<FixedImage, MovingImage, float> Registration;

    typedef itk::ImageFileReader<FixedImage> FixedImageReader;
    typedef itk::ImageFileReader<MovingImage> MovingImageReader;

    typedef itk::QuantileThresholdImageCalculator<FixedImage> FixedThresCalculator;
    typedef itk::QuantileThresholdImageCalculator<MovingImage>
    MovingThresCalculator;
    typedef itk::StatisticsImageFilter<FixedImage>     FixedStatsCalculator;
    typedef itk::StatisticsImageFilter<MovingImage>    MovingStatsCalculator;

    typedef itk::IntensityWindowingImageFilter<FixedImage,FixedImage> FixedWindowFilter;
    typedef itk::IntensityWindowingImageFilter<MovingImage,MovingImage> MovingWindowFilter;

    bool useCoincidenceW;
    bool useAdaptativeHistograms;
    double gaussian_deviation;
    bool output_images = false;
    unsigned int mres_levels;
    unsigned int mi_bins;
    unsigned int l1_samples, l2_samples;
    bool first_level_swap = false;  // TODO full support: let swap every level
    bool useRats = false;
    bool noLog = false;
    bool useVPH, useOldVPH = false;
    bool useAllRange = false;
    bool useQuantiles = false;
    bool reorient= false;
    bool useMasks=false;
    bool useHistMatching = false;
    bool forceOrient = false;

    QuantileOption fixedQuantiles = FIXED_QUANTILES;
    QuantileOption movingQuantiles = MOVING_QUANTILES;

    std::string fixedImageName;
    std::string movingImageName;
    std::string fixedMaskName;
    std::string outPrefix;
    unsigned int init_type, tf_type;

    bpo::options_description desc("Usage");
    desc.add_options()
                    ("help", "show help message")
                    ("fixed-image", bpo::value< std::string >(&fixedImageName)->required(), "fixed image file")
                    ("moving-image", bpo::value< std::string >(&movingImageName)->required(), "moving image file")
                    ("fixed-mask,F", bpo::value< std::string >(&fixedMaskName), "mask for fixed image file")
                    ("levels,l",bpo::value< unsigned int >(&mres_levels)->default_value(MRES_LEVELS), "number of multi-resolution levels" )
                    ("bins,b",bpo::value< unsigned int >(&mi_bins)->default_value(MI_BINS), "number of bins for histogram computations of MI metrics" )
                    ("l1-samples",bpo::value< unsigned int >(&l1_samples)->default_value(NUM_SAMPLES), "number of random samples for level 1" )
                    ("l2-samples",bpo::value< unsigned int >(&l2_samples)->default_value(NUM_SAMPLES), "number of random samples for level 2" )
                    ("output-prefix,o", bpo::value< std::string >(&outPrefix), "prefix for output files")
                    ("output-all", bpo::bool_switch(&output_images), "output intermediate images")
                    ("no-output-moving", bpo::bool_switch(), "prevent writing resampled moving image output")
                    ("first-level-swap", bpo::bool_switch(&first_level_swap)->default_value(false), "swap moving & fixed image on first level")
                    ("gaussian-deviation,g", bpo::value< double >(&gaussian_deviation)->default_value(GAUSSIAN_FILTER_DEVIATION), "Gaussian bluring filter deviation width (mm)")
                    ("init-type,I", bpo::value< unsigned int >(&init_type)->default_value(Registration::THelper::INIT_MOMENTS), "Initialization type (0=NONE, 1=GEOMETRY, 2=MOMENTS, 3=AFFINE, 4=PA, 5=CUSTOM_TRANSFORM)" )
                    ("transform-type,T", bpo::value< unsigned int >(&tf_type)->default_value(Registration::THelper::VERSOR_RIGID_3D), "Transform type (0=Euler3D, 1=CenteredEuler3D, 2=VersorRigid3D, 3=VersorLinear3D, 4=Affine7p, 5=Affine9p, 6=Affine15p, 7=bspline)" )
                    ("reorient", bpo::bool_switch(&reorient)->default_value(false), "reorient images if they have different orientations")
                    ("force-orientation", bpo::bool_switch(&forceOrient)->default_value(false), "force moving image orientation to the one of fixed image")
                    ("no-log", bpo::bool_switch(&noLog)->default_value(false), "do not log performance")
                    ("use-vph", bpo::bool_switch(&useVPH)->default_value(false), "Use VPHTk RM-SPECT brain co-registration preprocessing")
                    ("use-old-vph", bpo::bool_switch(&useOldVPH)->default_value(false), "Use old VPHTk RM-SPECT brain co-registration preprocessing")
                    ("use-all-range", bpo::bool_switch(&useAllRange)->default_value(false), "do not consider every pixel below the mean as background")
                    ("fixed-quantiles", bpo::value<QuantileOption>(&fixedQuantiles)->multitoken(), "use specified quantiles for pixel range definition of fixed image")
                    ("moving-quantiles", bpo::value<QuantileOption>(&movingQuantiles)->multitoken(), "use specified quantiles for pixel range definition of fixed image")
                    ("use-masks", bpo::bool_switch(&useMasks)->default_value(false), "use available masks for registration")
                    ("match-histogram", bpo::bool_switch(&useHistMatching)->default_value(false), "match moving image histogram to fixed image histogram")
                    ("use-RATS", bpo::bool_switch(&useRats)->default_value(false), "do not use Robust Automated Threshold Selection (RATS)")
                    ("use-adapt-histograms", bpo::bool_switch(&useAdaptativeHistograms)->default_value(false), "use adaptative histograms")
                    ("use-coincidence-weighting,w", bpo::bool_switch(&useCoincidenceW)->default_value(false), "use coincidence weighting");
    bpo::positional_options_description pdesc;
    pdesc.add("fixed-image", 1);
    pdesc.add("moving-image", 2);
    bpo::variables_map vmap;
    bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vmap);

    if (vmap.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    try {
        bpo::notify(vmap);
    } catch ( boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::required_option> > &err){
        std::cout << "Error: " << err.what() << std::endl;
        std::cout << desc << std::endl;
        return EXIT_FAILURE;
    }

    clock_t preProcessStart = clock();

    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    FixedImageReader::Pointer fixedImageReader = FixedImageReader::New();
    MovingImageReader::Pointer movingImageReader = MovingImageReader::New();

    fixedImageReader->SetFileName(fixedImageName);
    movingImageReader->SetFileName(movingImageName);

    FixedImage::ConstPointer fixedImage = fixedImageReader->GetOutput();
    MovingImage::ConstPointer movingImage = movingImageReader->GetOutput();


    FixedImage::ConstPointer fixedMask;
    if ( ! fixedMaskName.empty() ){
        FixedImageReader::Pointer fixedMaskReader = FixedImageReader::New();
        fixedMaskReader->SetFileName(fixedMaskName);

        try {
            fixedImageReader->Update();
            movingImageReader->Update();
        } catch (itk::ExceptionObject & err) {
            std::cout << "ExceptionObject caught !" << std::endl;
            std::cout << err << std::endl;
            return EXIT_FAILURE;
        }

        try {
            fixedMaskReader->Update();
        } catch (...) {

        }
        fixedMask = fixedMaskReader->GetOutput();
    }

    if ( reorient ){
        MovingImage::Pointer mtmp = movingImageReader->GetOutput();
        mtmp->SetDirection( fixedImage->GetDirection() );
        if (output_images)
            SaveImageToFile<MovingImage> ( MovingImage::ConstPointer(mtmp), outPrefix + "movingImageDirection.nii.gz");

        typedef itk::OrientImageFilter< MovingImage, MovingImage >                Orienter;
        Orienter::Pointer o_m = Orienter::New();
        o_m->SetInput( mtmp );
        // o_m->SetDesiredCoordinateDirection( fixedImage->GetDirection() );

        if(output_images){
            o_m->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP );
            o_m->Update();
            SaveImageToFile<MovingImage> (o_m->GetOutput(), outPrefix + "movingReOrientedImageRIP.nii.gz");

            o_m->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI );
            o_m->Update();
            SaveImageToFile<MovingImage> (o_m->GetOutput(), outPrefix + "movingReOrientedImagePRI.nii.gz");

            o_m->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR );
            o_m->Update();
            SaveImageToFile<MovingImage> (o_m->GetOutput(), outPrefix + "movingReOrientedImageIPR.nii.gz");
        }

        o_m->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI );
        o_m->Update();
        movingImage = o_m->GetOutput();

        if(output_images) SaveImageToFile<MovingImage> (o_m->GetOutput(), outPrefix + "movingReOrientedImagePLI.nii.gz");
    }else if ( forceOrient ) {
        MovingImage::Pointer mtmp = movingImageReader->GetOutput();
        mtmp->SetDirection( fixedImage->GetDirection() );
        movingImage = mtmp;
        if (output_images)
            SaveImageToFile<MovingImage> ( movingImage , outPrefix + "movingImageDirection.nii.gz");

    }


    if(output_images){
        SaveImageToFile<FixedImage> (fixedImage, outPrefix + "fixedImage.nii.gz");
        SaveImageToFile<MovingImage> (movingImage, outPrefix + "movingImage.nii.gz");
        //      SaveImageToFile<FixedImage> (FixedImage::ConstPointer(fixedMask2), outPrefix + "fixedOrientedMask.hdr");
        if (fixedMask.IsNotNull())
            SaveImageToFile<FixedImage> (fixedMask, outPrefix + "fixedMask.nii.gz");

        if (fixedMask.IsNotNull()) {
            typedef itk::ResampleImageFilter<FixedImage,FixedImage> Res;
            Res::Pointer r= Res::New();
            r->SetInput(fixedMask);
            r->SetReferenceImage(fixedImage);
            r->SetUseReferenceImage(true);
            r->Update();
            SaveImageToFile<FixedImage> (FixedImage::ConstPointer(r->GetOutput()), outPrefix + "fixedOrientedMask.nii.gz");
        }
    }

    // ----------------------------------------------------------------
    // FIN BLOQUE DE LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------

    std::ofstream logfile((outPrefix + "log.txt").c_str());


    // ----------------------------------------------------------------
    // BLOQUE DE PREPARACIÓN DEL REGISTRO
    // ----------------------------------------------------------------
    Registration::Pointer registration = Registration::New();
    registration->SetMultiResolutionLevels(mres_levels);
    registration->SetOutputStream(&logfile);
    registration->SetSaveRegImages(output_images);
    registration->SetOutputDirectory(outPrefix);
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);

    if (fixedMask.IsNotNull())
        registration->SetFixedImageMask(fixedMask);

    if (first_level_swap)   {
        registration->SetFixedMovingSwapForLevel(0);
    }

    FixedImage::PixelType  f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMean();
    FixedImage::PixelType  f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMax();
    MovingImage::PixelType m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMean();
    MovingImage::PixelType m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMax();

    if ( useAllRange ) {
        f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMin();
        m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMin();
    }

    if ( useRats ) {
        //registration->GetFixedPreprocessHelper(0)->ComputeRATSThreshold(2.0, 10.0);
        //registration->GetMovingPreprocessHelper(0)->ComputeRATSThreshold(2.0, 10.0);

        registration->GetFixedPreprocessHelper(0)->ApplyRATSMaskToInput(1);
        registration->GetMovingPreprocessHelper(0)->ApplyRATSMaskToInput(2);

        if( output_images ){
            SaveImageToFile<FixedImage> ( registration->GetFixedPreprocessHelper(0)->GetRATSThresholdedOutput(1),
                    outPrefix + "fixedRATSImage.nii.gz");

            SaveImageToFile<MovingImage> (registration->GetMovingPreprocessHelper(0)->GetRATSThresholdedOutput(2),
                    outPrefix + "movingRATSImage.nii.gz");
        }
    }


    useQuantiles = ( vmap.count("fixed-quantiles") || vmap.count("moving-quantiles") );

    if ( vmap.count("fixed-quantiles") || useVPH ) {
            registration->GetFixedPreprocessHelper(0)->ComputeQuantileThresholds(fixedQuantiles.lo_thres.get(), fixedQuantiles.hi_thres.get());
            f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMaxThreshold();
            f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMinThreshold();
    }

    if ( vmap.count("moving-quantiles") || useVPH ) {
            registration->GetMovingPreprocessHelper(0)->ComputeQuantileThresholds(movingQuantiles.lo_thres.get(), movingQuantiles.hi_thres.get());
            m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMaxThreshold();
            m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMinThreshold();
    }


    // Logs initialization

    if (!noLog) {
        for( unsigned int i = 0; i<mres_levels; i++ ) {
            std::stringstream logFileName;
            logFileName << outPrefix << "opt_log_" << i << ".txt";
            std::ofstream *logfileOpt = new std::ofstream( logFileName.str().c_str() );
            registration->GetOptimizerHelper(i)->SetOutputStream( logfileOpt );
        }
    }


    // Level 0 processing
    if( useOldVPH ) {

        registration->GetFixedPreprocessHelper(0)->ComputeQuantileThresholds(fixedQuantiles.lo_thres.get(), fixedQuantiles.hi_thres.get());
        registration->GetMovingPreprocessHelper(0)->ComputeQuantileThresholds(movingQuantiles.lo_thres.get(), movingQuantiles.hi_thres.get());
        f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMinThreshold();
        f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMaxThreshold();
        m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMinThreshold();
        m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMaxThreshold();



        if (useOldVPH || useRats) {
            //registration->GetFixedPreprocessHelper(0)->ComputeRATSThreshold(2.0, 10.0);
            //registration->GetMovingPreprocessHelper(0)->ComputeRATSThreshold(2.0, 10.0);


            registration->GetFixedPreprocessHelper(0)->ApplyRATSMask(1);
            registration->GetMovingPreprocessHelper(0)->ApplyRATSMask(2);

            f_MinTh = 1.0;
            m_MinTh = 1.0;

            registration->SetMovingImageMask( registration->GetMovingPreprocessHelper(0)->GetRATSThresholdedOutput(2) );
            //registration->SetInitializationUseMasks();


            if( output_images ){
                SaveImageToFile<FixedImage> ( registration->GetFixedPreprocessHelper(0)->GetRATSThresholdedOutput(1),
                        outPrefix + "fixedRATSImage.nii.gz");

                SaveImageToFile<MovingImage> (registration->GetMovingPreprocessHelper(0)->GetRATSThresholdedOutput(2),
                        outPrefix + "movingRATSImage.nii.gz");
            }
        }
        registration->GetMovingPreprocessHelper(0)->SetStripEmptyHighBins();
        registration->GetFixedPreprocessHelper(0)->SetStripEmptyHighBins();

        if( output_images ){
            registration->GetFixedPreprocessHelper(0)->SaveBinarizedOutput( outPrefix + "fixedQuantileBinarized.nii.gz");
            registration->GetMovingPreprocessHelper(0)->SaveBinarizedOutput( outPrefix + "movingQuantileBinarized.nii.gz");
        }
    }

    if ( useHistMatching ) {
        registration->GetFixedPreprocessHelper(0)->ComputeQuantileThresholds(fixedQuantiles.lo_thres.get(), fixedQuantiles.hi_thres.get());
        registration->GetFixedPreprocessHelper(0)->SetStripEmptyHighBins();
        f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMaxThreshold();

        registration->GetMovingPreprocessHelper(0)->SetHistogramMatchingImage( registration->GetFixedPreprocessHelper(0)->GetOutputImage() );

        f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMean();
        registration->GetFixedPreprocessHelper(0)->SetMinThreshold( f_MinTh );

        m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMean();
        registration->GetMovingPreprocessHelper(0)->SetMinThreshold(m_MinTh);

        m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMaxThreshold();

        if( output_images ){
            registration->GetFixedPreprocessHelper(0)->SaveOutput( outPrefix + "fixedPreprocessHighStrip.nii.gz");
            registration->GetMovingPreprocessHelper(0)->SaveOutput( outPrefix + "movingPreprocessMatched.nii.gz");
        }
    }

    if( !noLog ) {
        logfile << "* Levels=" << mres_levels << std::endl;
        logfile << "* SaveRegImages=" << output_images << std::endl;
        logfile << "* FirstLevelSwap=" << first_level_swap << std::endl;
        logfile << "* Smooth=" << gaussian_deviation << std::endl;
        logfile << "* UseOldVPH?=" << useOldVPH << std::endl;
        logfile << "* Use RATS?=" << useRats << std::endl;
        logfile << "* UseAdaptativeHistograms?=" << useAdaptativeHistograms << std::endl;
        logfile << "* UseCoincidenceWeighting?=" << useCoincidenceW << std::endl;
        logfile << "* UseQuantileThresholds?=" << useQuantiles << std::endl;
        if ( vmap.count("fixed-quantiles") )
        logfile << "     - Fixed:  (" << fixedQuantiles.lo_thres.get() << "," << fixedQuantiles.hi_thres.get() << ")" << std::endl;
        if ( vmap.count("moving-quantiles") )
        logfile << "     - Moving: (" << movingQuantiles.lo_thres.get() << "," << movingQuantiles.hi_thres.get() << ")" << std::endl;
        logfile << "* L1 Samples?=" << l1_samples << std::endl;
        logfile << "* L2 Samples?=" << l2_samples << std::endl;
        logfile << "* InitializationType=" << init_type << std::endl;
        logfile << "* Fixed Image Data:" << std::endl;
        logfile << "     - Image Origin: " << fixedImage->GetOrigin() << std::endl;
        logfile << "     - Image Direction: " << fixedImage->GetDirection() << std::endl;
        logfile << "     - Image Range: (" << registration->GetFixedPreprocessHelper(0)->GetMin() << "," << registration->GetFixedPreprocessHelper(0)->GetMax() << ")" << std::endl;
        logfile << "     - Registration Range: (" << f_MinTh << "," << f_MaxTh << ")."  << std::endl;
        logfile << "* Moving Image Data:" << std::endl;
        logfile << "     - Image Origin: " << movingImage->GetOrigin() << std::endl;
        logfile << "     - Image Direction: " << movingImage->GetDirection() << std::endl;
        logfile << "     - Image Range: (" << registration->GetMovingPreprocessHelper(0)->GetMin() << "," << registration->GetMovingPreprocessHelper(0)->GetMax() << ")" << std::endl;
        logfile << "     - Registration Range: (" << m_MinTh << "," << m_MaxTh << ")."  << std::endl;

        if (fixedMask.IsNotNull()) {
            logfile << "* Fixed Mask Data:" << std::endl;
            logfile << "     - Image Origin: " << fixedMask->GetOrigin() << std::endl;
            logfile << "     - Image Direction: " << fixedMask->GetDirection() << std::endl;
        }
    }

    registration->GetTransformHelper(0)->SetTransformType( Registration::THelper::TRANSFORM_TYPE (tf_type));
    registration->GetTransformHelper(0)->SetInitializationType(Registration::THelper::INITIALIZATION_TYPE (init_type) );
    registration->SetInitializationUseMasks();
    registration->SetMetricUseMasks(0,useMasks);
    registration->SetUseBinarizedInit(true);
    // registration->SetInitSize( 30 );

    registration->SmoothImages(0, gaussian_deviation * 4, gaussian_deviation * 4); // 2.0, 2.5
    registration->ResampleImages(0, 64, 64);
    registration->GetOptimizerHelper(0)->SetType(Registration::OHelper::SPSA);
    registration->GetMetricHelper(0)->SetMetricType(Registration::MHelper::NMI);

    unsigned int npar = registration->GetTransformHelper(0)->GetNumberOfParameters();
    Registration::OHelper::ScalesType scales( npar );

    if ( scales.GetSize() == 6 ) {
        scales.Fill(1.0);
        scales[0] = scales[1] = scales[2] = 1.0 / (1.0 * DEG); // 1.0 / (1.5 * DEG);
    } else {
        scales.Fill(2500);
        scales[0]=scales[1]=scales[2]=5000;
        scales[3]=scales[4]=scales[5]= 1;
    }

    //  scales = registration->GetAutoScales(fixedImage.GetPointer(), false);

    Registration::OHelper::ParametersType o(registration->GetOptimizerHelper(0)->GetNumberOfParameters());
    o[0] = SPSA_MAX_ITERATIONS;
    o[1] = true;
    o[2] = SPSA_a;
    o[3] = SPSA_A;
    o[4] = SPSA_C;
    o[5] = SPSA_MIN_ITERATIONS;
    o[8] = SPSA_NUM_PERTURBATIONS;
    o[9] = SPSA_INITIAL_STEP;
    registration->GetOptimizerHelper(0)->SetParameters(o);

    registration->GetOptimizerHelper(0)->SetScales(scales);

    Registration::MHelper::ParametersType p = registration->GetMetricHelper(0)->GetParameters();
    Registration::MHelper::NMIType::MeasurementVectorType lBound;
    Registration::MHelper::NMIType::MeasurementVectorType uBound;
    lBound.SetSize(2);
    uBound.SetSize(2);

    // MRI Bounds
    lBound[0] = f_MinTh; uBound[0] = f_MaxTh;
    // SPECT Bounds
    lBound[1] = m_MinTh; uBound[1] = m_MaxTh;

    if (first_level_swap){
        // MRI Bounds
        lBound[1] = f_MinTh; uBound[1] = f_MaxTh;
        // SPECT Bounds
        lBound[0] = m_MinTh; uBound[0] = m_MaxTh;
    }

    Registration::MHelper::MIType::HistogramType::SizeType histogramSize;
    histogramSize.SetSize(2);
    histogramSize.Fill( mi_bins );

    p[0] = histogramSize;
    p[1] = lBound;
    p[2] = uBound;
    p[3] = (double) lBound[0];
    p[5] = useAdaptativeHistograms;

    if( l2_samples > 0 ) {
        p[6]= l2_samples;
        p[8]=false;
    } else {
        p[8] = true;
    }

    registration->GetMetricHelper(0)->SetParameters(p);

    if ( mres_levels > 1 ) {
        // Level 1 processing
        registration->GetFixedPreprocessHelper(1)->SetMaxThreshold (f_MaxTh);
        registration->GetFixedPreprocessHelper(1)->SetMinThreshold (f_MinTh);
        registration->GetMovingPreprocessHelper(1)->SetMaxThreshold(m_MaxTh);
        registration->GetMovingPreprocessHelper(1)->SetMinThreshold(m_MinTh);
        registration->SetMetricUseMasks(1,useMasks);

        if ( useQuantiles || useVPH ) {
            registration->GetFixedPreprocessHelper(1)->SetStripEmptyHighBins( );
            registration->GetMovingPreprocessHelper(1)->SetStripEmptyHighBins( );
        }

        if ( useOldVPH ) {
            registration->GetFixedPreprocessHelper(1)->SetStripEmptyHighBins();
            registration->GetMovingPreprocessHelper(1)->SetStripEmptyHighBins();

            if (useOldVPH || useRats) {
                registration->GetFixedPreprocessHelper(1)->SetRATSThreshold(registration->GetFixedPreprocessHelper(0)->GetRATSThreshold());
                registration->GetFixedPreprocessHelper(1)->ApplyRATSMask(1);
                registration->GetMovingPreprocessHelper(1)->SetRATSThreshold(registration->GetMovingPreprocessHelper(0)->GetRATSThreshold());
                registration->GetMovingPreprocessHelper(1)->ApplyRATSMask(2);
            }
        }

        registration->SmoothImages(1, gaussian_deviation, gaussian_deviation * 2);

        if( output_images ){
            SaveImageToFile<Registration::FixedPHelper::OutputImageType> ( registration->GetFixedPreprocessHelper(1)->GetOutputImage(),
                    outPrefix + "fixedEqualizedImage.nii.gz");
            SaveImageToFile<Registration::MovingPHelper::OutputImageType> ( registration->GetMovingPreprocessHelper(1)->GetOutputImage(),
                    outPrefix + "movingEqualizedImage.nii.gz");
        }

        registration->GetOptimizerHelper(1)->SetType(Registration::OHelper::ROBBINS_MONRO);

        registration->GetTransformHelper(1)->SetTransformType(Registration::THelper::VERSOR_RIGID_3D);
        if ( npar > registration->GetTransformHelper(1)->GetNumberOfParameters() ) {
            // TODO throw exception
            registration->GetTransformHelper(1)->SetTransformType( registration->GetTransformHelper(0)->GetType() );
            logfile << "* Warning: changing level 1 transform type to match number of parameters" << std::endl;
        }

        registration->GetMetricHelper(1)->SetMetricType(Registration::MHelper::MATTES);

//      Registration::OHelper::ScalesType spsa_sc(6);
//      spsa_sc = registration->GetAutoScales(fixedImage.GetPointer(), true);

        Registration::OHelper::ParametersType o_mattes =
                registration->GetOptimizerHelper(1)->GetParameters();
        o_mattes[0] = GD_MAX_ITERATIONS;
        o_mattes[1] = false;
        o_mattes[2] = GD_A;
        o_mattes[3] = GD_a;
        o_mattes[4] = GD_ALPHA;
        registration->GetOptimizerHelper(1)->SetParameters(o_mattes);
        registration->GetOptimizerHelper(1)->SetScales(scales);

        Registration::MHelper::ParametersType p_mattes =
                registration->GetMetricHelper(1)->GetParameters();
        p_mattes[0] = static_cast<unsigned long> ( mi_bins );
        p_mattes[1] = l2_samples;
        p_mattes[5] = (double) f_MinTh;
        registration->GetMetricHelper(1)->SetParameters(p_mattes);
    }

    clock_t preProcessStop = clock();
    float pre_tot_t = (float) (((double)(preProcessStop - preProcessStart)) / CLOCKS_PER_SEC);
    int h, min, sec;
    h = (pre_tot_t/3600);
    min = (((int)pre_tot_t)%3600)/60;
    sec = (((int)pre_tot_t)%3600)%60;

    char pre_time[50];
    sprintf(pre_time, "\t* Pre-processing Total Time = %02d:%02d:%02d hours\n", h, min, sec);


    clock_t initTime = clock();

    int result = registration->StartRegistration();

    clock_t finishTime = clock();
    float tot_t = static_cast<float>( static_cast<double>(finishTime - initTime) / CLOCKS_PER_SEC);
    h = (tot_t/3600);
    min = (static_cast<int>(tot_t)%3600)/60;
    sec = (static_cast<int>(tot_t)%3600)%60;

    char reg_time[50];
    sprintf(reg_time, "\t* Registration Total Time = %02d:%02d:%02d hours\n", h, min, sec);


    logfile << pre_time << std::endl;
    logfile << reg_time << std::endl;

    registration->GetTransformHelper(mres_levels-1)->SaveTransformToFile( outPrefix + "tforms.txt" );

    if( vmap["no-output-moving"].as<bool>() == false)
        registration->GetTransformHelper(mres_levels-1)->ResampleImage< MovingImage, FixedImage >( MovingImage::ConstPointer (movingImage), FixedImage::ConstPointer (fixedImage), false, outPrefix + "movingRegistered.nii.gz" );

    return result;
}

void validate(boost::any& v, const std::vector<std::string>& values,
              QuantileOption*, int) {
    QuantileOption model;

    model.default_opt = boost::lexical_cast<bool> (values.at(0));

    if ( ! model.default_opt ) {
        if (values.size() >= 3){
            float hi_thres = boost::lexical_cast<float>(values.at(2));

            if ( hi_thres > 0.0 && hi_thres <= 1.0 ) {
                model.default_opt = false;
                model.hi_thres = hi_thres;
            } else {
                throw bpo::validation_error( bpo::validation_error::invalid_option_value, "Error specifying quantiles: Invalid upper quantile specified");
            }
        }

        if (values.size() >= 2) {
            float lo_thres = boost::lexical_cast<float>(values.at(1));

            if ( !model.default_opt && model.hi_thres <= lo_thres ) {
                throw bpo::validation_error( bpo::validation_error::invalid_option_value, "Error specifying quantiles: Lower quantile is greater than upper one");
            }

            if ( lo_thres >= 0.0 && lo_thres <= 1.0) {
                model.default_opt = false;
                model.lo_thres = lo_thres;
            } else {
                throw bpo::validation_error( bpo::validation_error::invalid_option_value, "Error specifying quantiles: Invalid lower quantile specified");
            }
        }
    }

    v = model;
}
