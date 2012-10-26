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

#ifndef MIREGISTRATION_TXX
#define MIREGISTRATION_TXX

#ifndef NUMBER_OF_THREADS
#define NUMBER_OF_THREADS 4
#endif

#include "MIRegistration.h"

template<class TFixedImage, class TMovingImage, class TScalarType>
MIRegistration<TFixedImage, TMovingImage, TScalarType>::MIRegistration() :
	m_UseMasksInitialization(false), m_OutputDirectory("./"), m_SaveRegImages(false), m_UseBinarizedInit(true), m_UseFixedOriginalInit(true), m_UseMovingOriginalInit(true), m_OutputStream(&std::cout) {
	// Initialize 1 level default configuration
	this->SetMultiResolutionLevels(1);
	m_InitFixedLowerBound = itk::NumericTraits<InputFixedImagePixelType>::min();
	m_InitMovingLowerBound = itk::NumericTraits<InputMovingImagePixelType>::min();
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetMultiResolutionLevels(unsigned int levels) {
	m_Levels = levels;
	m_InitSize = 0;

	// Initialize vectors:
	m_FixedPreprocessHelpers.resize(m_Levels);
	m_MovingPreprocessHelpers.resize(m_Levels);
	m_OptimizerHelpers.resize(m_Levels);
	m_MetricHelpers.resize(m_Levels);
	m_TransformHelpers.resize(m_Levels);
	m_RegistrationPointers.resize(m_Levels);
	m_FixedImage.resize(m_Levels);
	m_MovingImage.resize(m_Levels);
	m_SwapLevels.resize(m_Levels);
	m_UseMasksRegistrationLevel.resize(m_Levels);

	for (unsigned int i = 0; i < m_Levels; i++) {
		m_UseMasksRegistrationLevel[i] = false;
		m_SwapLevels[i] = false;

		if (m_FixedPreprocessHelpers[i].IsNull())
			m_FixedPreprocessHelpers[i] = FixedPHelper::New();
		if (m_MovingPreprocessHelpers[i].IsNull())
			m_MovingPreprocessHelpers[i] = MovingPHelper::New();
		if (m_OptimizerHelpers[i].IsNull())
			m_OptimizerHelpers[i] = OptimizerHelper::New();
		if (m_MetricHelpers[i].IsNull())
			m_MetricHelpers[i] = MHelper::New();
		if (m_TransformHelpers[i].IsNull())
			m_TransformHelpers[i] = THelper::New();
	}
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetInitialTransform(typename MIRegistration<TFixedImage, TMovingImage, TScalarType>::TransformPointer transform) {
	m_TransformHelpers[0]->SetInitialTransform(transform);
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetFixedImage(typename MIRegistration<TFixedImage, TMovingImage, TScalarType>::InputFixedImageConstPointer fixedImage) {
	m_InputFixedImage = fixedImage;
	for (unsigned int i = 0; i < m_Levels; i++)
		m_FixedPreprocessHelpers[i]->SetInputImage(fixedImage);
}
;

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetMovingImage(typename MIRegistration<TFixedImage, TMovingImage, TScalarType>::InputMovingImageConstPointer movingImage) {
	m_InputMovingImage = movingImage;
	for (unsigned int i = 0; i < m_Levels; i++)
		m_MovingPreprocessHelpers[i]->SetInputImage(movingImage);
}
;

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::NormalizeImages(unsigned int level, double minValue, double maxValue) {
	m_FixedPreprocessHelpers[level]->SetNormalization(minValue, maxValue);
	m_MovingPreprocessHelpers[level]->SetNormalization(minValue, maxValue);
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SmoothImages(unsigned int level, double smoothFixed, double smoothMoving) {
	if (smoothFixed > 0.0)
		m_FixedPreprocessHelpers[level]->SetGaussianSmoothFilter(smoothFixed);

	if (smoothMoving > 0.0)
		m_MovingPreprocessHelpers[level]->SetGaussianSmoothFilter(smoothMoving);
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::ResampleImages(unsigned int level, unsigned int fixedNOPD, unsigned int movingNOPD) {
	typename InternalFixedImageType::SizeType sizeFixed;
	sizeFixed.Fill(fixedNOPD);
	m_FixedPreprocessHelpers[level]->SetResample(sizeFixed);

	typename InternalFixedImageType::SizeType sizeMoving;
	sizeMoving.Fill(movingNOPD);
	m_MovingPreprocessHelpers[level]->SetResample(sizeMoving);
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetFixedMovingSwapForLevel(unsigned int level, bool value) {
	m_SwapLevels[level] = value;
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::Connect(unsigned int level) {
	if (m_RegistrationPointers[level].IsNull()) {
		m_RegistrationPointers[level] = RegistrationType::New();
		m_RegistrationPointers[level]->SetInterpolator(Interpolator::New());
	}

	if (m_SwapLevels[level]) {
		m_FixedImage[level] = m_MovingPreprocessHelpers[level]->GetOutputImage();
		m_MovingImage[level] = m_FixedPreprocessHelpers[level]->GetOutputImage();
	} else {
		m_MovingImage[level] = m_MovingPreprocessHelpers[level]->GetOutputImage();
		m_FixedImage[level] = m_FixedPreprocessHelpers[level]->GetOutputImage();
	}

	if (m_UseMasksRegistrationLevel[level]) {
		if (m_FixedMask.IsNotNull()) {
			typename InternalFixedCaster::Pointer caster = InternalFixedCaster::New();
			caster->SetInput(m_FixedMask);
			caster->Update();
			m_MetricHelpers[level]->SetFixedImageMask(caster->GetOutput());

		}
		if (m_MovingMask.IsNotNull()) {
			typename InternalMovingCaster::Pointer caster = InternalMovingCaster::New();
			caster->SetInput(m_MovingMask);
			caster->Update();
			m_MetricHelpers[level]->SetMovingImageMask(caster->GetOutput());
		}
	}

	m_RegistrationPointers[level]->SetFixedImage(m_FixedImage[level]);
	m_RegistrationPointers[level]->SetFixedImageRegion(m_FixedImage[level]->GetBufferedRegion());
	m_RegistrationPointers[level]->SetMovingImage(m_MovingImage[level]);
	m_RegistrationPointers[level]->SetNumberOfThreads(NUMBER_OF_THREADS);

	// Initializations
	if (level == 0) {
		InitImagesProcessing(m_UseBinarizedInit);
	} else {
		if (m_SwapLevels[level - 1]) {
			m_TransformHelpers[level]->SetInitialTransform(m_TransformHelpers[level - 1]->GetInverseTransform());
		} else {
			m_TransformHelpers[level]->SetInitialTransform(m_TransformHelpers[level - 1]->GetTransform());
		}
	}

	// Connect transformations
	m_TransformHelpers[level]->Connect(m_RegistrationPointers[level].GetPointer());
	// Connect metric
	m_MetricHelpers[level]->Connect(m_RegistrationPointers[level].GetPointer());
	// Connect optimizer
	m_OptimizerHelpers[level]->Connect(m_RegistrationPointers[level].GetPointer());

}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::InitImagesProcessing(bool binarize) {
	bool useFixedMask = m_UseMasksInitialization && m_FixedMask.IsNotNull();
	bool useMovingMask = m_UseMasksInitialization && m_MovingMask.IsNotNull();

	*m_OutputStream << "* Initialization: m_UseMasksInitialization?=" << m_UseMasksInitialization << std::endl;

	InternalFixedConstPointer fixed = m_FixedImage[0];
	if (useFixedMask) {
		*m_OutputStream << "* Initialization: Using Fixed Mask" << std::endl;
		typename InternalFixedCaster::Pointer cast = InternalFixedCaster::New();
		cast->SetInput(m_InputFixedImage);
		cast->Update();

		typename InternalFixedImageMasker::Pointer masker = InternalFixedImageMasker::New();
		masker->SetInput1(cast->GetOutput());
		masker->SetInput2(m_FixedMask);
		masker->Update();
		fixed = masker->GetOutput();
	}
	/*
	 typename InternalFixedCaster::Pointer cast = InternalFixedCaster::New();
	 cast->SetInput( m_InputFixedImage );
	 cast->Update();
	 fixed = cast->GetOutput();
	 */

	InternalMovingConstPointer moving = m_MovingImage[0];
	if (useMovingMask) {
		*m_OutputStream << "* Initialization: Using Moving Mask" << std::endl;
		typename InternalMovingCaster::Pointer cast = InternalMovingCaster::New();
		cast->SetInput(m_InputMovingImage);
		cast->Update();

		typename InternalMovingImageMasker::Pointer masker = InternalMovingImageMasker::New();
		masker->SetInput1(cast->GetOutput());
		masker->SetInput2(m_MovingMask);
		masker->Update();
		moving = masker->GetOutput();
	}

	/*
	 else if( m_UseMovingOriginalInit ){
	 *m_OutputStream << "* Initialization: Using Original Moving" << std::endl;
	 typename InternalMovingCaster::Pointer cast = InternalMovingCaster::New();
	 cast->SetInput( m_InputMovingImage );
	 cast->Update();
	 moving = cast->GetOutput();
	 }*/

	if (m_SaveRegImages) {
		SaveImageToFile<InternalFixedImageType> (fixed, m_OutputDirectory + "fixedFullSizeInit.nii.gz");
		SaveImageToFile<InternalMovingImageType> (moving, m_OutputDirectory + "movingFullSizeInit.nii.gz");
	}

	if (m_InitSize > 0) {
		// Initialize at a lower scale
		typename InternalFixedImageType::SizeType f_size = fixed->GetLargestPossibleRegion().GetSize();
		typename InternalMovingImageType::SizeType m_size = moving->GetLargestPossibleRegion().GetSize();
		unsigned int f_factors[3] = { 1, 1, 1 };
		unsigned int m_factors[3] = { 1, 1, 1 };
		bool f_resize = false;
		bool m_resize = false;

		for (unsigned int i = 0; i < FixedImageDimension; i++) {
			f_factors[i] = std::max((int) vcl_floor( (float) (f_size[i] / m_InitSize)), 1);
			m_factors[i] = std::max((int) vcl_floor( (float) (m_size[i] / m_InitSize)), 1);
			if (f_factors[i] > 1.0)
				f_resize = true;
			if (m_factors[i] > 1.0)
				m_resize = true;
		}

		if (f_resize) {
			typedef itk::ShrinkImageFilter<InternalFixedImageType, InternalFixedImageType> FixedShrink;
			typename FixedShrink::Pointer s = FixedShrink::New();
			s->SetInput(fixed);
			s->SetShrinkFactors(f_factors);
			s->Update();
			fixed = s->GetOutput();
		}

		if (m_resize) {
			typedef itk::ShrinkImageFilter<InternalMovingImageType, InternalMovingImageType> MovingShrink;
			typename MovingShrink::Pointer s = MovingShrink::New();
			s->SetInput(moving);
			s->SetShrinkFactors(m_factors);
			s->Update();
			moving = s->GetOutput();
		}

		if (m_SaveRegImages) {
			SaveImageToFile<InternalFixedImageType> (fixed, m_OutputDirectory + "fixedResizedInit.nii.gz");
			SaveImageToFile<InternalMovingImageType> (moving, m_OutputDirectory + "movingResizedInit.nii.gz");
		}
	}
	/*
	 else {
	 typedef itk::ResampleImageFilter< InternalFixedImageType, InternalFixedImageType> FixedResampler;
	 typedef itk::ResampleImageFilter< InternalMovingImageType, InternalMovingImageType> MovingResampler;

	 typename InternalFixedImageType::SizeType hresSize;
	 hresSize.

	 }
	 */

	typedef itk::StatisticsImageFilter<InternalFixedImageType> FixedCalc;
	typename FixedCalc::Pointer f_c = FixedCalc::New();
	f_c->SetInput(fixed);
	f_c->Update();

	if (m_InitFixedLowerBound <= itk::NumericTraits<InputFixedImagePixelType>::min()) {
		if (!useFixedMask) {
			m_InitFixedLowerBound = f_c->GetMean();
		} else
			m_InitFixedLowerBound += 0.1;
	}

	if (binarize) {
		typedef itk::BinaryThresholdImageFilter<InternalFixedImageType, InternalFixedImageType> FixedBinaryFilter;
		typename FixedBinaryFilter::Pointer f_b = FixedBinaryFilter::New();
		f_b->SetInput(fixed);
		f_b->SetLowerThreshold(m_InitFixedLowerBound);
		f_b->SetInsideValue(1);
		f_b->SetOutsideValue(0);
		f_b->Update();
		fixed = f_b->GetOutput();
	} else {
		typedef itk::IntensityWindowingImageFilter<InternalFixedImageType, InternalFixedImageType> FixedBinaryFilter;
		typename FixedBinaryFilter::Pointer f_b = FixedBinaryFilter::New();
		f_b->SetInput(fixed);
		f_b->SetWindowMinimum(m_InitFixedLowerBound);
		f_b->SetWindowMaximum(f_c->GetMaximum());
		f_b->SetOutputMinimum(0);
		f_b->SetOutputMaximum(300);
		f_b->Update();
		fixed = f_b->GetOutput();
	}

	typedef itk::StatisticsImageFilter<InternalMovingImageType> MovingCalc;
	typename MovingCalc::Pointer m_c = MovingCalc::New();
	m_c->SetInput(moving);
	m_c->Update();

	if (m_InitMovingLowerBound <= itk::NumericTraits<InputMovingImagePixelType>::min()) {
		if (!useMovingMask) {
			if (useFixedMask && binarize) {
				typename FixedCalc::Pointer f_m_c = FixedCalc::New();
				f_m_c->SetInput(fixed);
				f_m_c->Update();
				double percentage = f_m_c->GetMean();
				typename InternalFixedImageType::SpacingType sp = fixed->GetSpacing();
				double volume = percentage * fixed->GetLargestPossibleRegion().GetNumberOfPixels() * sp[0] * sp[1] * sp[2];
				typename InternalMovingImageType::SpacingType sp_m = moving->GetSpacing();
				double movingVolume = moving->GetLargestPossibleRegion().GetNumberOfPixels() * sp_m[0] * sp_m[1] * sp_m[2];

				double m_percentage = volume / movingVolume;

				typedef typename itk::Statistics::ScalarImageToHistogramGenerator<InternalMovingImageType> HistogramGenerator;
				typedef typename HistogramGenerator::HistogramType Histogram;
				typename HistogramGenerator::Pointer histogramGenerator = HistogramGenerator::New();
				histogramGenerator->SetInput(moving);
				histogramGenerator->SetNumberOfBins(256);
				histogramGenerator->SetMarginalScale(10.0);
				histogramGenerator->Compute();
				const Histogram * histogram = histogramGenerator->GetOutput();
				m_InitMovingLowerBound = histogram->Quantile(0, (1 - m_percentage));
			} else {
				m_InitMovingLowerBound = m_c->GetMean();
			}
		}
	}

	if (binarize) {
		typedef itk::BinaryThresholdImageFilter<InternalMovingImageType, InternalMovingImageType> MovingBinaryFilter;
		typename MovingBinaryFilter::Pointer m_b = MovingBinaryFilter::New();
		m_b->SetInput(moving);
		m_b->SetLowerThreshold(m_InitMovingLowerBound);
		m_b->SetInsideValue(1);
		m_b->SetOutsideValue(0);
		m_b->Update();
		moving = m_b->GetOutput();
	} else {
		typedef itk::IntensityWindowingImageFilter<InternalMovingImageType, InternalMovingImageType> MovingBinaryFilter;
		typename MovingBinaryFilter::Pointer m_b = MovingBinaryFilter::New();
		m_b->SetInput(moving);
		m_b->SetWindowMinimum(m_InitMovingLowerBound);
		m_b->SetWindowMaximum(m_c->GetMaximum());
		m_b->SetOutputMinimum(0);
		m_b->SetOutputMaximum(300);
		m_b->Update();
		moving = m_b->GetOutput();
	}

	if (m_SaveRegImages) {
		SaveImageToFile<InternalFixedImageType> (fixed, m_OutputDirectory + "fixedInit.nii.gz");
		SaveImageToFile<InternalMovingImageType> (moving, m_OutputDirectory + "movingInit.nii.gz");
	}

	m_TransformHelpers[0]->SetInitializationImages(fixed, moving, m_SwapLevels[0]);
}
;

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SaveResampledImage() {
	m_TransformHelpers[0]->ResampleImage<InputFixedImageType, InputMovingImageType> (m_InputMovingImage, m_InputFixedImage, m_SwapLevels[0], "res_0_pre.nii.gz");

}

template<class TFixedImage, class TMovingImage, class TScalarType>
int MIRegistration<TFixedImage, TMovingImage, TScalarType>::StartRegistration() {
	for (unsigned int i = 0; i < m_Levels; i++) {
		this->Connect(i);

		*m_OutputStream << "Level [" << i << "] Pre Helpers: ----------------" << std::endl;
		m_TransformHelpers[i]->Print(*m_OutputStream);
		m_OptimizerHelpers[i]->Print(*m_OutputStream);
		m_MetricHelpers[i]->Print(*m_OutputStream);

		std::stringstream name;
		name << m_OutputDirectory << "in_level_" << i;

		if (m_SaveRegImages) {
			m_TransformHelpers[i]->ResampleImage<InputMovingImageType, InputFixedImageType> (m_InputMovingImage, m_InputFixedImage, m_SwapLevels[i], name.str() + "_pre.nii.gz");
			SaveImageToFile<InternalFixedImageType> (m_FixedImage[i], name.str() + "_regFixedImage.nii.gz");
			SaveImageToFile<InternalMovingImageType> (m_MovingImage[i], (name.str() + "_regMovingImage.nii.gz"));
		}

		try {
			m_RegistrationPointers[i]->StartRegistration();
		} catch (itk::ExceptionObject & err) {
			*m_OutputStream << "ExceptionObject caught !" << std::endl;
			*m_OutputStream << err << std::endl;
			m_OptimizerHelpers[i]->Print(*m_OutputStream);
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::SetOutputStream(std::ostream* os) {
	m_OutputStream = os;
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void MIRegistration<TFixedImage, TMovingImage, TScalarType>::PrintSelf(std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	for (unsigned int i = 0; i < m_Levels; i++) {
		os << indent << "Components for level " << i << "/" << m_Levels << ":" << std::endl;
		m_FixedImage[i]->Print(os, indent);
		m_MovingImage[i]->Print(os, indent);
		m_OptimizerHelpers[i]->Print(os, indent);
		m_TransformHelpers[i]->Print(os, indent);
		m_MetricHelpers[i]->Print(os, indent);
	}
}

template<class TFixedImage, class TMovingImage, class TScalarType>
template<class TImage>
typename MIRegistration<TFixedImage, TMovingImage, TScalarType>::ScalesType MIRegistration<TFixedImage, TMovingImage, TScalarType>::GetAutoScales(const TImage* image, bool normalize) const {
	double center_index[3];
	typedef typename TImage::PointType Point;
	typedef itk::ImageRegionConstIterator<TImage> Iterator;
	typedef typename THelper::VersorRigidType Transform;

	ScalesType scales(6);
	scales.Fill(1.0);

	typename Transform::Pointer inc_t = Transform::New();
	typename Transform::Pointer dec_t = Transform::New();

	inc_t->SetIdentity();
	dec_t->SetIdentity();

	typename TImage::SizeType imSize = image->GetLargestPossibleRegion().GetSize();
	typename TImage::SpacingType spacing = image->GetSpacing();
	Point origin = image->GetOrigin();

	center_index[0] = (imSize[0] - 1) * 0.5;
	center_index[1] = (imSize[1] - 1) * 0.5;
	center_index[2] = (imSize[2] - 1) * 0.5;

	Point center_point;
	center_point[0] = center_index[0] * spacing[0] + origin[0];
	center_point[1] = center_index[1] * spacing[1] + origin[1];
	center_point[2] = center_index[2] * spacing[2] + origin[2];

	inc_t->SetCenter(center_point);
	dec_t->SetCenter(center_point);

	double volume = (imSize[0] * spacing[0]) * (imSize[1] * spacing[1]) * (imSize[2] * spacing[2]);
	double inc_v = spacing[0] * spacing[1] * spacing[2];
	double v = 0;
	double max = 0;

	typename THelper::GenericType::ParametersType initialParameters = inc_t->GetParameters();

	for (unsigned int index = 0; index < scales.size(); index++) {
		double norm = 0;

		typename THelper::GenericType::ParametersType inc_param = initialParameters;
		typename THelper::GenericType::ParametersType dec_param = initialParameters;

		double adjustUnit = (index <= 2) ? DEG : 1.0;
		double increment = 1.0 * adjustUnit;
		inc_param[index] += (increment * 0.5);
		dec_param[index] -= (increment * 0.5);

		inc_t->SetParameters(inc_param);
		dec_t->SetParameters(dec_param);

		Iterator it(image, image->GetLargestPossibleRegion());
		it.GoToBegin();
		while (!it.IsAtEnd()) {

			Point p;
			image->TransformIndexToPhysicalPoint(it.GetIndex(), p);

			Point inc_p = inc_t->TransformPoint(p);
			Point dec_p = dec_t->TransformPoint(p);

			Point inc;
			inc.Fill(0.0);
			inc[0] = (inc_p[0] - dec_p[0]) / increment;
			inc[1] = (inc_p[1] - dec_p[1]) / increment;
			inc[2] = (inc_p[2] - dec_p[2]) / increment;

			double distance = (pow((double) inc[0], (int) 2u) + pow((double) inc[1], (int) 2u) + pow((double) inc[2], (int) 2u));

			if (distance > 0)
				norm += distance * inc_v;

			v += inc_v;
			++it;
		}

		scales[index] = sqrt(norm / volume);

		if (scales[index] > max) {
			max = scales[index];
		}

	}

	if (normalize) {
		for (unsigned int i = 0; i < scales.size(); i++)
			scales[i] /= max;
	}

	return scales;
}

#endif // MIREGISTRATION_TXX
