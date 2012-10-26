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



typedef  float  PixelRange;
const unsigned int Dimension = 3;
typedef itk::OrientedImage< PixelRange, Dimension >          FixedImage;
typedef itk::OrientedImage< PixelRange, Dimension >          MovingImage;

const double GAUSSIAN_FILTER_DEVIATION =        1.0;
const double IMAGES_MAX_VALUE          =      300.0;

#ifndef USE_ADAPTATIVE_BINING
const double SPSA_A                      =    100.0; // 500.0; // 10.0
const double SPSA_C                      =      0.1; // 0.8;   // 0.08
const double SPSA_a                      =    10000; // 10000; // 2000
const unsigned int SPSA_MAX_ITERATIONS   =     4500; // 5000;  // 2500
#else
const double SPSA_A                      =    450.0;
const double SPSA_C                      =      0.1;
const double SPSA_a                      =     1000;
const unsigned int SPSA_MAX_ITERATIONS   =     4500;
#endif

const double GD_INITIAL_LEARNING_RATE    =    25e-2;
const unsigned int GD_MAX_ITERATIONS     =     1000;
const double  GD_ALPHA                   =     0.95;
const double  GD_a                       =     30.0;
const double  GD_A                       =    100.0;


const bool OUTPUT_IMAGES                 =    true;



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
            if (lastBinUnderMean == i-1 ) {  binsUnderMean++;   }
            else                          {  binsUnderMean = 1; }
            
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
	if( argc < 3 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFileName  movingImageFileName ";
		std::cerr << " [fixedMaskFileName] ";
		std::cerr << " [outputDirectoryName=./]";
		std::cerr << "SetUseCoincidenceWeighting";
		std::cerr <<  std::endl;
		return EXIT_FAILURE;
	}



	typedef MIRegistration< FixedImage, MovingImage,float >     Registration;

	// Lectores. Se crean dos readers, porque la imagen fija y móvil
	// pueden ser diferentes
	typedef itk::ImageFileReader< FixedImage  > FixedImageReader;
	typedef itk::ImageFileReader< MovingImage > MovingImageReader;

// 	// Escritores
// 	typedef itk::ImageFileWriter< FixedImage >  FixedWriter;
// 	typedef itk::ImageFileWriter< MovingImage >  MovingWriter;



	// ----------------------------------------------------------------
	// LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------
	FixedImageReader::Pointer  fixedImageReader  = FixedImageReader::New();
	MovingImageReader::Pointer movingImageReader = MovingImageReader::New();

	std::string fixedImageName  = argv[1];
	std::string movingImageName = argv[2];
	std::string fixedMaskName   = argv[3];
	std::string outputDirName   = "./";
	
	if( argc >= 5 )
	{
	  outputDirName = argv[4];
	}
	
	
	fixedImageReader->SetFileName(  fixedImageName );
	movingImageReader->SetFileName( movingImageName );
	
	FixedImageReader::Pointer fixedMaskReader = FixedImageReader::New();
	fixedMaskReader->SetFileName( fixedMaskName );
	

    FixedImage::ConstPointer fixedMask;
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

    try {
        fixedMaskReader->Update();
        fixedMask = fixedMaskReader->GetOutput();
    }
    catch (...)
    {

    }

 	FixedImage::ConstPointer fixedImage = fixedImageReader->GetOutput();
	MovingImage::ConstPointer movingImage = movingImageReader->GetOutput();
	
// 	typedef itk::OrientImageFilter< FixedImage, FixedImage >  Orienter;
// 	Orienter::Pointer orienter = Orienter::New();
// 	orienter->SetInput( fixedImageReader->GetOutput() );
// 	orienter->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS);
// 	orienter->Update();
// 	fixedImage = orienter->GetOutput();
/*
	typedef itk::OrientImageFilter< MovingImage, MovingImage >  OrienterMoving;
	OrienterMoving::Pointer orienterMoving = OrienterMoving::New();
	orienterMoving->SetInput( movingImageReader->GetOutput() );
	orienterMoving->SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS);
	orienterMoving->Update();
	movingImage = orienterMoving->GetOutput();*/
	
	// ----------------------------------------------------------------
	// FIN BLOQUE DE LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------

	
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
	registration->SetMultiResolutionLevels(2);
	registration->SetFixedMovingSwapForLevel(0);
	
	registration->GetOptimizerHelper(0)->SetType( Registration::OHelper::SPSA );
	registration->GetTransformHelper(0)->SetTransformType( Registration::THelper::VERSOR_RIGID_3D );
	registration->SetInitializationUseMasks();
        registration->GetMetricHelper(0)->SetMetricType( Registration::MHelper::NMI );
	
	registration->GetOptimizerHelper(1)->SetType( Registration::OHelper::ROBBINS_MONRO );
	registration->GetTransformHelper(1)->SetTransformType( Registration::THelper::VERSOR_RIGID_3D );
	registration->GetMetricHelper(   1)->SetMetricType( Registration::MHelper::MATTES );
	
	Registration::OHelper::ScalesType scales ( 6 );

	scales = registration->GetAutoScales( fixedImage.GetPointer(), false );
	Registration::OHelper::ParametersType o (registration->GetOptimizerHelper(0)->GetNumberOfParameters());
	o[0] = SPSA_MAX_ITERATIONS; o[1] = true; o[2]= SPSA_a; o[3] = SPSA_A; o[4] = SPSA_C;
	registration->GetOptimizerHelper(0)->SetParameters(o);
	registration->GetOptimizerHelper(0)->SetScales( scales );

	registration->SetFixedImage( fixedImage ); //, fixedThPercent );
	if ( fixedMask.IsNotNull() ) registration->SetFixedImageMask( fixedMask );

	registration->SetMovingImage( movingImage ); //, movingThPercent );
	
	// Level 0 processing
	registration->NormalizeImages( 0, 0.0, IMAGES_MAX_VALUE );
	registration->SmoothImages( 0, GAUSSIAN_FILTER_DEVIATION * 6, GAUSSIAN_FILTER_DEVIATION * 8 );
	registration->ResampleImages( 0, 64, 64 );
	
	// Level 1 processing
	registration->NormalizeImages( 1, 0.0, IMAGES_MAX_VALUE );
	registration->SmoothImages( 1, GAUSSIAN_FILTER_DEVIATION * 1.5, GAUSSIAN_FILTER_DEVIATION * 2 );
	
	Registration::MHelper::ParametersType p = registration->GetMetricHelper(0)->GetParameters();
	
#ifndef USE_ADAPTATIVE_BINING
	Registration::MHelper::NMIType::MeasurementVectorType lBound;
	Registration::MHelper::NMIType::MeasurementVectorType uBound;
	lBound.SetSize( 2 );
	uBound.SetSize( 2 );
	
	// MRI Bounds
	lBound[1] = 15.0;
	uBound[1] = IMAGES_MAX_VALUE * 1.001;

	// SPECT Bounds
	lBound[0] = 50.0;
	uBound[0] = IMAGES_MAX_VALUE * 1.001;
	
	p[1] = lBound;
	p[2] = uBound;
	p[3] = (double) lBound[0];
#else
	p[4] = true;
	p[5] = true;
	p[6] = 40000u;
#endif

	registration->GetMetricHelper(0)->SetParameters(p);
		
	Registration::OHelper::ParametersType o_mattes = registration->GetOptimizerHelper(1)->GetParameters();
	o_mattes[0] = GD_MAX_ITERATIONS;
	o_mattes[1] = false;
	o_mattes[2] = GD_A;
	o_mattes[3] = GD_a;
	o_mattes[4] = GD_ALPHA;
	registration->GetOptimizerHelper(1)->SetParameters( o_mattes );
	registration->GetOptimizerHelper(1)->SetScales( scales );
	
	Registration::MHelper::ParametersType p_mattes = registration->GetMetricHelper(1)->GetParameters();
	p_mattes[5] = 15.0;
	registration->GetMetricHelper(1)->SetParameters(p_mattes);
	
	registration->SetSaveRegImages( OUTPUT_IMAGES );
	registration->SetOutputDirectory( outputDirName );

	
	SaveImageToFile<FixedImage>( fixedImage, outputDirName + "fixedImage.nii" );
	
	std::ofstream logfile( (outputDirName + "log.txt").c_str() );
 	registration->SetOutputStream( &logfile );

	return registration->StartRegistration();

}
