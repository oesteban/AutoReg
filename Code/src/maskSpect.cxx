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


// #define USE_ADAPTATIVE_BINING
#include "MIRegistration.h"


typedef short PixelRange;
const unsigned int Dimension = 3;
typedef itk::Image<PixelRange, Dimension> FixedImage;
typedef itk::Image<PixelRange, Dimension> MovingImage;

const unsigned int MATTES_SPATIAL_SAMPLES = 7000;

const double GAUSSIAN_FILTER_DEVIATION = 40.0;
const double FIXED_IMAGE_BG_PERCENTAGE  = 0.70;
const double FIXED_IMAGE_MAX_PERCENTAGE = 0.999;
const double MOVING_IMAGE_BG_PERCENTAGE = 0.70;
const double MOVING_IMAGE_MAX_PERCENTAGE= 0.99;

const unsigned int GD_MAX_ITERATIONS =  2000;
const double GD_ALPHA = 0.95;
const double GD_a = 200.0;
const double GD_A = 400;

const bool OUTPUT_IMAGES = true;


int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " templateImage targetImage maskTemplate outputFileName ";
		std::cerr << " [outputDirectoryName=./]";
		std::cerr << "SetUseCoincidenceWeighting";
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	typedef MIRegistration<FixedImage, MovingImage, float> Registration;

	// Lectores. Se crean dos readers, porque la imagen fija y móvil
	// pueden ser diferentes
	typedef itk::ImageFileReader<FixedImage> FixedImageReader;
	typedef itk::ImageFileReader<MovingImage> MovingImageReader;

	typedef itk::QuantileThresholdImageCalculator<FixedImage> FixedThresCalculator;
	typedef itk::QuantileThresholdImageCalculator<MovingImage>
				MovingThresCalculator;
	typedef itk::StatisticsImageFilter<FixedImage>     FixedStatsCalculator;
	typedef itk::StatisticsImageFilter<MovingImage>	   MovingStatsCalculator;

	typedef itk::IntensityWindowingImageFilter<FixedImage,FixedImage> FixedWindowFilter;
	typedef itk::IntensityWindowingImageFilter<MovingImage,MovingImage> MovingWindowFilter;
	// ----------------------------------------------------------------
	// LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------
	FixedImageReader::Pointer fixedImageReader = FixedImageReader::New();
	MovingImageReader::Pointer movingImageReader = MovingImageReader::New();

	std::string fixedImageName = argv[1];
	std::string movingImageName = argv[2];
	std::string maskImageName = argv[3];
	std::string outputFileName = argv[4];
	std::string outputDirName = "./";

	if (argc >= 6) {
		outputDirName = argv[argc-1];
	}

	std::ofstream logfile((outputDirName + "log.txt").c_str());

	fixedImageReader->SetFileName(fixedImageName);
	movingImageReader->SetFileName(movingImageName);

	try {
		fixedImageReader->Update();
		movingImageReader->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}


	FixedImage::ConstPointer fixedImage = fixedImageReader->GetOutput();
	MovingImage::ConstPointer movingImage = movingImageReader->GetOutput();

	// ----------------------------------------------------------------
	// FIN BLOQUE DE LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------


	// ----------------------------------------------------------------
	// BLOQUE DE PREPARACIÓN DEL REGISTRO
	// ----------------------------------------------------------------
	Registration::Pointer registration = Registration::New();
	registration->SetMultiResolutionLevels(1);

	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);

	FixedImage::PixelType  f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMean();
	FixedImage::PixelType  f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMax();
	MovingImage::PixelType m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMean();
	MovingImage::PixelType m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMax();

/*
	registration->GetFixedPreprocessHelper(0)->ComputeQuantileThresholds(FIXED_IMAGE_BG_PERCENTAGE,FIXED_IMAGE_MAX_PERCENTAGE);
	registration->GetMovingPreprocessHelper(0)->ComputeQuantileThresholds(MOVING_IMAGE_BG_PERCENTAGE,MOVING_IMAGE_MAX_PERCENTAGE);

	FixedImage::PixelType f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMinThreshold();
	//FixedImage::PixelType f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMaxThreshold();
	MovingImage::PixelType m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMinThreshold();
	//MovingImage::PixelType m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMaxThreshold();

	//registration->SetUseBinarizedInit(false);
	//registration->SetUseFixedOriginalInit(false);
	//registration->SetUseMovingOriginalInit(false);
	registration->SetInitFixedLowerBound(f_MinTh);
	registration->SetInitMovingLowerBound(m_MinTh);*/

	registration->GetFixedPreprocessHelper(0)->ApplyRATSMaskToInput(2);
	registration->GetMovingPreprocessHelper(0)->ApplyRATSMaskToInput(2);


	// Level 0 processing
	registration->GetFixedPreprocessHelper(0)->SetStripEmptyHighBins();
	registration->GetMovingPreprocessHelper(0)->SetStripEmptyHighBins();

	if( OUTPUT_IMAGES ){
		SaveImageToFile<FixedImage> ( registration->GetFixedPreprocessHelper(0)->GetRATSThresholdedOutput(2),
				outputDirName + "fixedRATSImage.nii.gz");

		SaveImageToFile<MovingImage> (registration->GetMovingPreprocessHelper(0)->GetRATSThresholdedOutput(2),
				outputDirName + "movingRATSImage.nii.gz");
	}

	f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMean();
	f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMax();
	m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMean();
	m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMax();
	registration->SetUseBinarizedInit(true);


	if( OUTPUT_IMAGES ){
		SaveImageToFile<Registration::FixedPHelper::OutputImageType> (
				registration->GetFixedPreprocessHelper(0)->GetOutputImage(),
				outputDirName + "fixedEqualizedImage.nii");
		SaveImageToFile<Registration::MovingPHelper::OutputImageType> (
				registration->GetMovingPreprocessHelper(0)->GetOutputImage(),
				outputDirName + "movingEqualizedImage.nii");
	}

	registration->SmoothImages(0, 0.0, GAUSSIAN_FILTER_DEVIATION);

	registration->GetOptimizerHelper(0)->SetType(Registration::OHelper::ROBBINS_MONRO);
	registration->GetTransformHelper(0)->SetTransformType(Registration::THelper::AFFINE_9P);
	registration->GetTransformHelper(0)->SetInitializationType(Registration::THelper::INIT_AFFINE);
	registration->GetMetricHelper(0)->SetMetricType(Registration::MHelper::MATTES);

	Registration::OHelper::ScalesType scales(9);
	scales.Fill(2500);
	scales[0]=scales[1]=scales[2]=5000;
	scales[3]=scales[4]=scales[5]= 1;
	//scales = registration->GetAutoScales(fixedImage.GetPointer(), false);
	registration->GetOptimizerHelper(0)->SetScales(scales);

	Registration::OHelper::ParametersType o(
			registration->GetOptimizerHelper(0)->GetNumberOfParameters());

	o[0] = GD_MAX_ITERATIONS;
	o[1] = false;
	o[2] = GD_A;
	o[3] = GD_a;
	o[4] = GD_ALPHA;
	registration->GetOptimizerHelper(0)->SetParameters(o);


	std::ofstream opt_0_logfile((outputDirName + "optimizer_log_0.txt").c_str());
	registration->GetOptimizerHelper(0)->SetOutputStream( &opt_0_logfile );

	Registration::MHelper::ParametersType p = registration->GetMetricHelper(0)->GetParameters();
	p[1] = MATTES_SPATIAL_SAMPLES;
	p[5] = (double) f_MinTh;
	registration->GetMetricHelper(0)->SetParameters(p);

	registration->SetSaveRegImages(OUTPUT_IMAGES);
	registration->SetOutputDirectory(outputDirName);

	if( OUTPUT_IMAGES )
	  SaveImageToFile<FixedImage> (fixedImage, outputDirName + "fixedImage.nii");

	registration->SetOutputStream(&logfile);

	int result = registration->StartRegistration();

	FixedImageReader::Pointer maskReader = FixedImageReader::New();
	maskReader->SetFileName(maskImageName);
	try {
		maskReader->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	registration->GetTransformHelper(0)->ResampleImage< FixedImage, MovingImage >( FixedImage::ConstPointer (maskReader->GetOutput()), movingImage, true, outputDirName + outputFileName );

	if( OUTPUT_IMAGES )
		registration->GetTransformHelper(0)->ResampleImage< FixedImage, MovingImage >(fixedImage,movingImage, true, outputDirName + "templateRegistered.nii.gz" );

	return result;


}
