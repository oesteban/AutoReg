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

//#define USE_ADAPTATIVE_BINING
#include "MIRegistration.h"

typedef  short  FixedPixelRange;
typedef  float  MovingPixelRange;
const unsigned int Dimension = 3;
typedef itk::Image< FixedPixelRange, Dimension >          FixedImage;
typedef itk::Image< MovingPixelRange, Dimension >         MovingImage;

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
        std::cerr << " fixedImageFile  movingImageFile [mattes]";
        std::cerr <<  std::endl;
        return EXIT_FAILURE;
    }


    typedef MIRegistration< FixedImage, MovingImage,float >     Registration;
    typedef itk::ImageFileReader< FixedImage  > FixedImageReader;
    typedef itk::ImageFileReader< MovingImage > MovingImageReader;
    typedef itk::ImageFileWriter< FixedImage >  FixedWriter;
    typedef itk::ImageFileWriter< MovingImage >  MovingWriter;

    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    FixedImageReader::Pointer  fixedImageReader  = FixedImageReader::New();
    MovingImageReader::Pointer movingImageReader = MovingImageReader::New();

    fixedImageReader->SetFileName(  argv[1] );
    movingImageReader->SetFileName( argv[2] );
    
    std::string metric_name( argv[3] );
    
    bool is_nmi = metric_name.compare("mattes")!=0;

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

    // ----------------------------------------------------------------
    // FIN BLOQUE DE LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    
    typedef Registration::THelper::EulerType      EulerTf;
    EulerTf::Pointer tf = EulerTf::New();
    tf->SetIdentity();
    
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

    registration->SetFixedImage( fixedImage ); //, fixedThPercent );
    registration->SetMovingImage( movingImage.GetPointer() ); // ,movingThPercent );

    Registration::MHelper::ParametersType p;

  if(!is_nmi)
  {
    registration->GetMetricHelper(0)->SetMetricType( Registration::MHelper::MATTES );
    p = registration->GetMetricHelper(0)->GetParameters();
    p[5] = 15.0;
  }
  else
  {
    registration->GetMetricHelper(0)->SetMetricType( Registration::MHelper::NMI );

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

    p = registration->GetMetricHelper(0)->GetParameters();
    p[1] = lBound;
    p[2] = uBound;
    p[3] = 15.0;
#endif

    p[7] = 8u;
  }

    registration->GetMetricHelper(0)->SetParameters(p);
    
    
#ifndef USE_ADAPTATIVE_BINING
    registration->NormalizeImages( 0, 0.0, IMAGES_MAX_VALUE );
#endif
    registration->SmoothImages( 0, GAUSSIAN_FILTER_DEVIATION * 1.5, GAUSSIAN_FILTER_DEVIATION * 2 );
    
    Registration::OHelper::ParametersType o (registration->GetOptimizerHelper(0)->GetNumberOfParameters());
    Registration::OHelper::ExhaustiveType::StepsType steps ( 6 );
    steps.Fill( 10 );
    o[0] = steps;

    Registration::OHelper::ScalesType scales ( 6 );
    double rotationRadians = ( M_PI / 360 );
    scales[0] = rotationRadians;
    scales[1] = rotationRadians;
    scales[2] = rotationRadians;
    scales[3] = 1.0;
    scales[4] = 1.0;
    scales[5] = 1.0;
    registration->GetOptimizerHelper(0)->SetParameters(o);
    registration->GetOptimizerHelper(0)->SetScales( scales );

    return registration->StartRegistration();
}
