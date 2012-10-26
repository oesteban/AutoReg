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



#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkGIBUB3DTransform.h>

//#define USE_ADAPTATIVE_BINING
#include "MIRegistration.h"

typedef  float  PixelRange;
const unsigned int Dimension = 3;
typedef itk::Image< PixelRange, Dimension >          FixedImage;
typedef itk::Image< PixelRange, Dimension >          MovingImage;

const double GAUSSIAN_FILTER_DEVIATION =        1.0;
const double IMAGES_MAX_VALUE          =      300.0;
const double EXPECTED_OFFSET             =     20.0; // 30 // in mm
const double EXPECTED_ROTATION_MAGNITUDE =     10.0; // in degrees (7.0)
const double GD_EXPECTED_OFFSET          =      4.0; // in mm
const double GD_EXPECTED_ROTATION_MAGNITUDE =   1.7; // in degrees (0.7)
const double SPSA_A                      =      100;
const double SPSA_C                      =     0.08; // 0.1
const double SPSA_a                      =     40.0; // 16! // 0.5 0.16




template <class TInputImage >
void ImageHistogramEqualization( typename TInputImage::Pointer image, typename TInputImage::Pointer out, float highThres )
{
    const unsigned int HISTOGRAM_NUMBER_OF_BINS = 256;

    typedef typename TInputImage::PixelType                                InputPixel;
    typedef itk::MinimumMaximumImageCalculator< TInputImage >              InputImageCalculator;
    typedef itk::IntensityWindowingImageFilter< TInputImage, TInputImage > InputRescale;
    typedef itk::ThresholdImageFilter< TInputImage >                       InputThreshold;

    typename InputImageCalculator::Pointer calculator = InputImageCalculator::New();
    calculator->SetImage(image);
    calculator->Compute();
    InputPixel img_min = calculator->GetMinimum();
    InputPixel img_max = calculator->GetMaximum();

    typedef typename itk::Statistics::ScalarImageToHistogramGenerator<TInputImage>   HistogramGenerator;
    typedef typename HistogramGenerator::HistogramType  Histogram;
    typename HistogramGenerator::Pointer histogramGenerator = HistogramGenerator::New();
    histogramGenerator->SetInput( image );
    histogramGenerator->SetNumberOfBins( HISTOGRAM_NUMBER_OF_BINS );
    histogramGenerator->SetMarginalScale( 10.0 );
    histogramGenerator->Compute();
    const Histogram * histogram = histogramGenerator->GetOutput();

    typename Histogram::BinMaxContainerType maxs = histogram->GetMaxs();
    typename Histogram::BinMinContainerType mins = histogram->GetMins();

    typename Histogram::BinMaxVectorType binDimensionMaximums = histogram->GetDimensionMaxs( 0 );
    const unsigned int histogramSize = histogram->Size();

    typename Histogram::AbsoluteFrequencyType maxFreq = 0;
    unsigned int maxBin = 0;
    int minBin = -1;
    const unsigned int minAllowedFreq = highThres * (histogram->GetTotalFrequency() / HISTOGRAM_NUMBER_OF_BINS);

    unsigned int binsUnderMean = 0;
    unsigned int lastBinUnderMean = 0;

    for (unsigned int i = 0; i<histogramSize; i++ )
    {
        typename Histogram::AbsoluteFrequencyType freq = histogram->GetFrequency( i , 0 );

        if (i>0 && freq>maxFreq)
        {
            maxFreq = freq;
            maxBin = i;
        }

        if (i> 0.25* HISTOGRAM_NUMBER_OF_BINS && freq < minAllowedFreq && minBin == -1)
        {
            if (lastBinUnderMean == i-1 ) {
                binsUnderMean++;
            }
            else                          {
                binsUnderMean = 1;
            }

            lastBinUnderMean = i;

            if ( binsUnderMean >= 5 ) minBin = i-binsUnderMean+1;
        }
    }

    InputPixel binMax = histogram->GetBinMax( 0, maxBin );
    InputPixel binMin = histogram->GetBinMin( 0, maxBin );
    InputPixel intensityMax = (binMax+binMin)* 0.5;

    InputPixel window_min = img_min;
    InputPixel window_max = img_max;

    if ( maxBin < (HISTOGRAM_NUMBER_OF_BINS / 4.0) )
        window_min = intensityMax;

    if (minBin > -1 ) {
        binMax = histogram->GetBinMax( 0, minBin );
        binMin = histogram->GetBinMin( 0, minBin );
        window_max = 0.5*( binMax + binMin );
    }

    typename InputThreshold::Pointer thres = InputThreshold::New();
    thres->SetInput( image );
    thres->SetOutsideValue( 0 );
    thres->ThresholdOutside( window_min, window_max );
    thres->Update();

    typename InputRescale::Pointer rescaler = InputRescale::New();
    rescaler->SetOutputMinimum(img_min);
    rescaler->SetOutputMaximum(img_max);
    rescaler->SetInput( thres->GetOutput() );
    rescaler->SetWindowMinimum(window_min);
    rescaler->SetWindowMaximum(window_max);
    rescaler->Update();

    out = rescaler->GetOutput();
}

void FocusDetFixedInit( FixedImage::Pointer fixed )
{
    FixedImage::SizeType size = fixed->GetLargestPossibleRegion().GetSize();

    if ( size[0] > 256 || size[1] > 256 || size[2] > 256 )
    {
        FixedImage::SizeType newSize = size;
        FixedImage::SpacingType spacing = fixed->GetSpacing();
        FixedImage::SpacingType newSpacing = spacing;
        for ( unsigned int i = 0; i < FixedImage::ImageDimension; i++ )
        {
            if ( newSize[i] > 256 ) {
                newSize[i] = 256;
                newSpacing[i] = (size[i]*spacing[i])/newSize[i];
            }
        }

        // TODO: Sustituir por un shrink filter
        typedef itk::ResampleImageFilter< FixedImage, FixedImage > Resampler;
        Resampler::Pointer fRes = Resampler::New();
        fRes->SetInput( fixed );
        fRes->SetOutputDirection( fixed->GetDirection() );
        fRes->SetOutputOrigin( fixed->GetOrigin() );
        fRes->SetOutputSpacing( newSpacing );
        fRes->SetSize( newSize );
        fRes->Update();
        fixed = fRes->GetOutput();
    }

    ImageHistogramEqualization< FixedImage >(fixed, fixed, 0.001);
    // End Fixed Image Initialization ------------------------------------------------------
};

void FocusDetMovingInit( MovingImage::Pointer moving ) {
    // Moving Image Initialization ------------------------------------------------------------
    ImageHistogramEqualization< MovingImage > ( moving, moving, 0.005 );
    // End Moving Image Initialization ------------------------------------------------------------
};

int main( int argc, char *argv[] )
{
    if ( argc < 2 )
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0];
        std::cerr << " fixedImageFile  movingImageFile ";
        std::cerr << "outputImagefile ";
        std::cerr << "NumberOfHistogramBinsForWritingTheMutualInformationHistogramMetric";
        std::cerr << "SetUseCoincidenceWeighting";
        std::cerr <<  std::endl;
        return EXIT_FAILURE;
    }



    typedef MIRegistration< FixedImage, MovingImage,float >     Registration;
    typedef itk::GIBUB3DTransform< double >                     Transform;

    // Lectores. Se crean dos readers, porque la imagen fija y móvil
    // pueden ser diferentes
    typedef itk::ImageFileReader< FixedImage  > FixedImageReader;
    typedef itk::ImageFileReader< MovingImage > MovingImageReader;

    // Escritores
    typedef itk::ImageFileWriter< FixedImage >  FixedWriter;
    typedef itk::ImageFileWriter< MovingImage >  MovingWriter;

    double angles[3];
    double translations[3];

    unsigned int case_id = atoi(argv[1]);
    std::string data_dir = "/home/oesteban/workspace/oesteban/GIBUBTransformTest/testingData/";
    std::stringstream rm_dir;   rm_dir << data_dir << "RM" << std::setw(2) << std::setfill('0') <<  case_id << "/";
    
    std::string tf_file_name = data_dir + "gold-tforms.csv";
    std::ifstream gold_tf_file( tf_file_name.c_str() );

    if ( gold_tf_file.is_open() )
    {
        std::string s;
        gold_tf_file.seekg (0, std::ios::beg); // go to the first line

        for (unsigned int i=1; i< case_id; i++) // loop 'till the desired line
            std::getline(gold_tf_file, s);

        std::getline( gold_tf_file, s, ',' );
        if ( static_cast<unsigned int>(atoi(s.c_str()))!=case_id )
        {
            std::cout << "Bad line" << std::endl;
            return EXIT_FAILURE;
        }

        std::getline( gold_tf_file, s, ',' );
        translations[0] = atof( s.c_str() );
        std::getline( gold_tf_file, s, ',' );
        translations[1] = atof( s.c_str() );
        std::getline( gold_tf_file, s, ',' );
        translations[2] = atof( s.c_str() );
        std::getline( gold_tf_file, s, ',' );
        angles[0] = atof( s.c_str() );
        std::getline( gold_tf_file, s, ',' );
        angles[1] = atof( s.c_str() );
        std::getline( gold_tf_file, s, ',' );
        angles[2] = atof( s.c_str() );
    }
    else
    {
        std::cout << "Gold Transforms file not found" << std::endl;
        return EXIT_FAILURE;
    }
    

    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    FixedImageReader::Pointer  fixedImageReader  = FixedImageReader::New();
    MovingImageReader::Pointer movingImageReader = MovingImageReader::New();

    std::string fixedImageName = rm_dir.str() + "xxRM3dT1x" + std::string( argv[1] ) + ".hdr";
    std::string movingImageName = rm_dir.str() + "BAS.trasz.centro.RM." + std::string( argv[1] ) + ".hdr";
    fixedImageReader->SetFileName(  fixedImageName.c_str() );
    movingImageReader->SetFileName( movingImageName.c_str() );

    try {
        fixedImageReader->Update();
        movingImageReader->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
    }

    FixedImage::ConstPointer fixedImage = fixedImageReader->GetOutput();
    MovingImage::Pointer movingImage = movingImageReader->GetOutput();


    FixedImage::SizeType      s_fixed  = fixedImage->GetLargestPossibleRegion().GetSize();
    FixedImage::IndexType     last_i_fixed;
    last_i_fixed[0] = s_fixed[0] - 1;
    last_i_fixed[1] = s_fixed[1] - 1;
    last_i_fixed[2] = s_fixed[2] - 1;

    FixedImage::PointType     o_fixed  = fixedImage->GetOrigin();
    FixedImage::SpacingType   sp_fixed = fixedImage->GetSpacing();

    MovingImage::SizeType      s_moving  = movingImage->GetLargestPossibleRegion().GetSize();
    MovingImage::IndexType last_i_moving;
    last_i_moving[0] = s_moving[0] - 1;
    last_i_moving[1] = s_moving[1] - 1;
    last_i_moving[2] = s_moving[2] - 1;

    MovingImage::PointType     o_moving  = movingImage->GetOrigin();
    MovingImage::SpacingType   sp_moving = movingImage->GetSpacing();
    
    MovingImage::PointType center_moving;
    double center_index_moving[3];
    center_index_moving[0] = last_i_moving[0] * 0.5;
    center_index_moving[1] = last_i_moving[1] * 0.5;
    center_index_moving[2] = last_i_moving[2] * 0.5;
    movingImage->TransformContinuousIndexToPhysicalPoint<double>( center_index_moving, center_moving );
    
    FixedImage::PointType center_fixed;
    double center_index_fixed[3];
    center_index_fixed[0] = (s_fixed[0]-1) * 0.5;
    center_index_fixed[1] = (s_fixed[1]-1) * 0.5;
    center_index_fixed[2] = (s_fixed[2]-1) * 0.5;
    fixedImage->TransformContinuousIndexToPhysicalPoint<double>( center_index_fixed , center_fixed );

    Registration::THelper::EulerType::InputVnlVectorType centers_vector;
    centers_vector[0] = (center_moving[0] - center_fixed[0]);
    centers_vector[1] = (center_moving[1] - center_fixed[1]);
    centers_vector[2] = (center_moving[2] - center_fixed[2]);

    MovingImage::PointType new_o_moving;
    new_o_moving[0] = o_moving[0] - centers_vector[0];
    new_o_moving[1] = o_moving[1] - centers_vector[1];
    new_o_moving[2] = o_moving[2] - centers_vector[2];
    
    movingImage->SetOrigin( new_o_moving );
    
    
    
    // ----------------------------------------------------------------
    // FIN BLOQUE DE LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    
    // Create and generate a Gold Registration transform.
    Transform::Pointer tf_gold     = Transform::New();
    tf_gold->SetIdentity();
    tf_gold->SetImageSpacing( sp_moving );
    tf_gold->SetReferenceSpacing( sp_fixed );
    tf_gold->SetGibUbCenter( center_index_moving );
    tf_gold->SetGibUbRotation( angles[0] * DEG , angles[1] * DEG ,angles[2] * DEG);
    tf_gold->SetGibUbTranslation( translations );
    
    typedef Registration::THelper::EulerType      EulerTf;
    EulerTf::Pointer tf = EulerTf::New();
    
    tf->SetCenter( tf_gold->GetCenter() );
    tf->SetRotation( tf_gold->GetAngleX(), tf_gold->GetAngleY(), tf_gold->GetAngleZ() );
    tf->SetTranslation( tf_gold->GetTranslation() );

    // ----------------------------------------------------------------
    // BLOQUE DE PREPARACIÓN DEL REGISTRO
    // ----------------------------------------------------------------
    Registration::Pointer   registration  = Registration::New();

    // Prepare Reg Image
    double fixedThPercent = 100.0;
    double movingThPercent = 100.0;
#ifndef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
    fixedThPercent = 35.0;
    movingThPercent = 7.0;
#endif
    registration->SetInitialTransform( tf );
    registration->SetMultiResolutionLevels(1);
    registration->GetOptimizerHelper(0)->SetType( Registration::OHelper::EXHAUSTIVE );
    registration->GetTransformHelper(0)->SetTransformType( Registration::THelper::EULER_3D );
    registration->GetMetricHelper(0)->SetMetricType( Registration::MHelper::NMI );



    Registration::OHelper::ParametersType o (registration->GetOptimizerHelper(0)->GetNumberOfParameters());
    Registration::OHelper::ExhaustiveType::StepsType steps ( 6 );
    steps.Fill( 12 );//50
    steps[0]=steps[1]=steps[2] = 40;
    o[0] = steps;
//     o[1] = 0.0175;


    Registration::OHelper::ScalesType scales ( 6 );
//    double rotationRadians = (2 * M_PI * EXPECTED_ROTATION_MAGNITUDE / 360 );
//    scales[0] = rotationRadians / EXPECTED_OFFSET;
//    scales[1] = rotationRadians / EXPECTED_OFFSET;
//    scales[2] = rotationRadians / EXPECTED_OFFSET;
    scales[0] = scales[1] = scales[2] = 0.02625;
    scales[3] = scales[4] = scales[5] = 4.0;

    registration->GetOptimizerHelper(0)->SetParameters(o);
    registration->GetOptimizerHelper(0)->SetScales( scales );


    registration->SetFixedImage( fixedImage ); // , GAUSSIAN_FILTER_DEVIATION * 3, fixedThPercent );
    registration->SetMovingImage( movingImage.GetPointer()); // , GAUSSIAN_FILTER_DEVIATION * 3, movingThPercent );

#ifndef USE_ADAPTATIVE_BINING
//    registration->NormalizeImages( 0, 0.0, IMAGES_MAX_VALUE );
#endif

    registration->ResampleImages( 0 , 50, 64 );

    Registration::MHelper::ParametersType p = registration->GetMetricHelper(0)->GetParameters();
    
#ifndef USE_ADAPTATIVE_BINING
    Registration::MHelper::NMIType::MeasurementVectorType lBound;
    Registration::MHelper::NMIType::MeasurementVectorType uBound;
    lBound.SetSize( 2 );
    uBound.SetSize( 2 );

    typedef itk::MinimumMaximumImageCalculator< FixedImage > FixedCalculator;
    FixedCalculator::Pointer c_f = FixedCalculator::New();
    c_f->SetImage ( fixedImage );
    c_f->Compute();
    // MRI Bounds
    lBound[1] = 15;
    uBound[1] = c_f->GetMaximum();

    typedef itk::MinimumMaximumImageCalculator< MovingImage > MovingCalculator;
    MovingCalculator::Pointer c_m = MovingCalculator::New();
    c_m->SetImage ( movingImage );
    c_m->Compute();
    // SPECT Bounds
    lBound[0] = 15;
    uBound[0] = c_m->GetMaximum();

    p[1] = lBound;
    p[2] = uBound;
    p[3] = 15.0;
#endif

    p[7] = 8u;

    registration->GetMetricHelper(0)->SetParameters(p);


    return registration->StartRegistration();
}
