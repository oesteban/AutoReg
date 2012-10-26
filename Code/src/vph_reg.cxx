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

#define USE_ADAPTATIVE_BINING
#include "MIRegistration.h"

typedef float PixelRange;
const unsigned int Dimension = 3;
typedef itk::Image<PixelRange, Dimension> FixedImage;
typedef itk::Image<PixelRange, Dimension> MovingImage;

const double GAUSSIAN_FILTER_DEVIATION = 1.0;
const double FIXED_IMAGE_BG_PERCENTAGE  = 0.65; // 0.55 - 0.75
const double FIXED_IMAGE_MAX_PERCENTAGE = 0.97; // 0.97
const double MOVING_IMAGE_BG_PERCENTAGE = 0.85; // 0.85
const double MOVING_IMAGE_MAX_PERCENTAGE= 0.999;

const unsigned int SPSA_MIN_ITERATIONS    = 100;
const unsigned int SPSA_NUM_PERTURBATIONS =   2;
const double       SPSA_INITIAL_STEP      = 1.0;

#ifndef USE_ADAPTATIVE_BINING
const double SPSA_A =                    50.0;
const double SPSA_C =                     0.4;
const double SPSA_a =                   1.0e4;
const unsigned int SPSA_MAX_ITERATIONS =  500;
#else
const double SPSA_A = 450.0;
const double SPSA_C = 0.1;
const double SPSA_a = 1000;
const unsigned int SPSA_MAX_ITERATIONS = 4500;
#endif

const double GD_INITIAL_LEARNING_RATE =  1000;
const unsigned int GD_MAX_ITERATIONS =  1000;
const double GD_ALPHA = 0.95;
const double GD_a = 30.0;
const double GD_A = 100;

const bool FIRST_LEVEL_SWAP = false;
const bool OUTPUT_IMAGES = false;


int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFileName  movingImageFileName ";
		std::cerr << " [fixedMaskFileName] ";
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
	std::string fixedMaskName = argv[3];
	std::string outputDirName = "./";

	if (argc >= 5) {
		outputDirName = argv[4];
	}

	std::ofstream logfile((outputDirName + "log.txt").c_str());

	fixedImageReader->SetFileName(fixedImageName);
	movingImageReader->SetFileName(movingImageName);

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

	FixedImage::ConstPointer fixedImage = fixedImageReader->GetOutput();
	MovingImage::ConstPointer movingImage = movingImageReader->GetOutput();
	FixedImage::ConstPointer fixedMask = fixedMaskReader->GetOutput();

	if(OUTPUT_IMAGES){
	    SaveImageToFile<FixedImage> (fixedImage, outputDirName + "fixedImage.nii.gz");
        SaveImageToFile<MovingImage> (movingImage, outputDirName + "movingImage.nii.gz");
  //      SaveImageToFile<FixedImage> (FixedImage::ConstPointer(fixedMask2), outputDirName + "fixedOrientedMask.hdr");
        SaveImageToFile<FixedImage> (fixedMask, outputDirName + "fixedMask.nii.gz");

		std::cout << "* Fixed Image Data:" << std::endl;
		std::cout << "     - Image Origin: " << fixedImage->GetOrigin() << std::endl;
		std::cout << "     - Image Direction: " << fixedImage->GetDirection() << std::endl;
		std::cout << "* Moving Image Data:" << std::endl;
		std::cout << "     - Image Origin: " << movingImage->GetOrigin() << std::endl;
		std::cout << "     - Image Direction: " << movingImage->GetDirection() << std::endl;
		std::cout << "* Fixed Mask Data:" << std::endl;
		std::cout << "     - Image Origin: " << fixedMask->GetOrigin() << std::endl;
		std::cout << "     - Image Direction: " << fixedMask->GetDirection() << std::endl;

		typedef itk::ResampleImageFilter<FixedImage,FixedImage> Res;
		Res::Pointer r= Res::New();
		r->SetInput(fixedMask);
		r->SetReferenceImage(fixedImage);
		r->SetUseReferenceImage(true);
		r->Update();
		SaveImageToFile<FixedImage> (FixedImage::ConstPointer(r->GetOutput()), outputDirName + "fixedOrientedMask.nii.gz");

	}

	// ----------------------------------------------------------------
	// FIN BLOQUE DE LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------


	// ----------------------------------------------------------------
	// BLOQUE DE PREPARACIÓN DEL REGISTRO
	// ----------------------------------------------------------------
	Registration::Pointer registration = Registration::New();
	registration->SetMultiResolutionLevels(2);

	if (FIRST_LEVEL_SWAP)	registration->SetFixedMovingSwapForLevel(0);

	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImage);

	if (fixedMask.IsNotNull())
		registration->SetFixedImageMask(fixedMask);

	registration->GetFixedPreprocessHelper(0)->ComputeQuantileThresholds(FIXED_IMAGE_BG_PERCENTAGE,FIXED_IMAGE_MAX_PERCENTAGE);
	registration->GetMovingPreprocessHelper(0)->ComputeQuantileThresholds(MOVING_IMAGE_BG_PERCENTAGE,MOVING_IMAGE_MAX_PERCENTAGE);
	registration->GetFixedPreprocessHelper(1)->ComputeQuantileThresholds(FIXED_IMAGE_BG_PERCENTAGE,FIXED_IMAGE_MAX_PERCENTAGE);
	registration->GetMovingPreprocessHelper(1)->ComputeQuantileThresholds(MOVING_IMAGE_BG_PERCENTAGE,MOVING_IMAGE_MAX_PERCENTAGE);

	FixedImage::PixelType f_MinTh = registration->GetFixedPreprocessHelper(0)->GetMinThreshold();
	FixedImage::PixelType f_MaxTh = registration->GetFixedPreprocessHelper(0)->GetMaxThreshold();
	MovingImage::PixelType m_MinTh = registration->GetMovingPreprocessHelper(0)->GetMinThreshold();
	MovingImage::PixelType m_MaxTh = registration->GetMovingPreprocessHelper(0)->GetMaxThreshold();
	
	

	if( OUTPUT_IMAGES ){
		registration->GetFixedPreprocessHelper(0)->SaveBinarizedOutput( outputDirName + "fixedQuantileBinarized.nii.gz");
		registration->GetMovingPreprocessHelper(0)->SaveBinarizedOutput( outputDirName + "movingQuantileBinarized.nii.gz");
	}

	// Level 0 processing
	registration->GetFixedPreprocessHelper(1)->SetStripEmptyHighBins();
	registration->GetMovingPreprocessHelper(0)->SetStripEmptyHighBins();
	registration->GetMovingPreprocessHelper(1)->SetStripEmptyHighBins();

	if( OUTPUT_IMAGES ){
		SaveImageToFile<Registration::FixedPHelper::OutputImageType> (
				registration->GetFixedPreprocessHelper(1)->GetOutputImage(),
				outputDirName + "fixedEqualizedImage.nii.gz");
		SaveImageToFile<Registration::MovingPHelper::OutputImageType> (
				registration->GetMovingPreprocessHelper(1)->GetOutputImage(),
				outputDirName + "movingEqualizedImage.nii.gz");
	}

	registration->SmoothImages(0, GAUSSIAN_FILTER_DEVIATION * 6, GAUSSIAN_FILTER_DEVIATION * 8);
	registration->ResampleImages(0, 64, 64);

	// Level 1 processing
	registration->SmoothImages(1, GAUSSIAN_FILTER_DEVIATION * 1.5,	GAUSSIAN_FILTER_DEVIATION * 2);


	registration->GetOptimizerHelper(0)->SetType(Registration::OHelper::SPSA);
	registration->GetTransformHelper(0)->SetTransformType(
			Registration::THelper::VERSOR_RIGID_3D);
	registration->GetTransformHelper(0)->SetInitializationType(Registration::THelper::INIT_MOMENTS);
	registration->SetInitializationUseMasks();
	registration->GetMetricHelper(0)->SetMetricType(Registration::MHelper::NMI);

	registration->GetOptimizerHelper(1)->SetType(
			Registration::OHelper::ROBBINS_MONRO);
	registration->GetTransformHelper(1)->SetTransformType(
			Registration::THelper::VERSOR_RIGID_3D);
	registration->GetMetricHelper(1)->SetMetricType(
			Registration::MHelper::MATTES);

	Registration::OHelper::ScalesType scales(6);

	scales = registration->GetAutoScales(fixedImage.GetPointer(), false);
	Registration::OHelper::ParametersType o(
			registration->GetOptimizerHelper(0)->GetNumberOfParameters());

	o[0] = SPSA_MAX_ITERATIONS;
	o[1] = true;
	o[2] = SPSA_a;
	o[3] = SPSA_A;
	o[4] = SPSA_C;
	o[5] = SPSA_MIN_ITERATIONS;
	o[8] = SPSA_NUM_PERTURBATIONS;
	o[9] = SPSA_INITIAL_STEP;
	registration->GetOptimizerHelper(0)->SetParameters(o);

	Registration::OHelper::ScalesType spsa_sc(6);
	spsa_sc = registration->GetAutoScales(fixedImage.GetPointer(), true);
	registration->GetOptimizerHelper(0)->SetScales(scales);

	std::ofstream opt_0_logfile((outputDirName + "optimizer_log_0.txt").c_str());
	std::ofstream opt_1_logfile((outputDirName + "optimizer_log_1.txt").c_str());
	registration->GetOptimizerHelper(0)->SetOutputStream( &opt_0_logfile );
	registration->GetOptimizerHelper(1)->SetOutputStream( &opt_1_logfile );


	Registration::MHelper::ParametersType p = registration->GetMetricHelper(0)->GetParameters();
	Registration::MHelper::NMIType::MeasurementVectorType lBound;
	Registration::MHelper::NMIType::MeasurementVectorType uBound;
	lBound.SetSize(2);
	uBound.SetSize(2);

	if (FIRST_LEVEL_SWAP){
		// MRI Bounds
		lBound[1] = f_MinTh; uBound[1] = f_MaxTh;
		// SPECT Bounds
		lBound[0] = m_MinTh; uBound[0] = m_MaxTh;
	}
	else {
		// MRI Bounds
		lBound[0] = f_MinTh; uBound[0] = f_MaxTh;
		// SPECT Bounds
		lBound[1] = m_MinTh; uBound[1] = m_MaxTh;
	}
	p[1] = lBound;
	p[2] = uBound;
	p[3] = (double) lBound[0];

	registration->GetMetricHelper(0)->SetParameters(p);

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
	p_mattes[5] = (double) f_MinTh;
	registration->GetMetricHelper(1)->SetParameters(p_mattes);

	registration->SetSaveRegImages(OUTPUT_IMAGES);
	registration->SetOutputDirectory(outputDirName);

	registration->SetOutputStream(&logfile);

	int result = registration->StartRegistration();

	registration->GetTransformHelper(1)->SaveTransformToFile( outputDirName + "tforms.txt" );
	registration->GetTransformHelper(1)->ResampleImage< MovingImage, FixedImage >( MovingImage::ConstPointer (movingImage), FixedImage::ConstPointer (fixedImage), false, outputDirName + "movingRegistered.nii.gz" );

	return result;
}
