/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkHistogramImageToImageMetric.txx,v $
 Language:  C++
 Date:      $Date: 2009-05-05 17:47:30 $
 Version:   $Revision: 1.30 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/

#ifndef __itkHistogramImageToImageMetric_txx
#define __itkHistogramImageToImageMetric_txx

#include "itkHistogramImageToImageMetric.h"

#include "itkArray.h"
#include "itkNumericTraits.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk {
template<class TFixedImage, class TMovingImage>
HistogramImageToImageMetric<TFixedImage, TMovingImage>::HistogramImageToImageMetric() {
	itkDebugMacro("Constructor");

	m_HistogramSize.Fill(256);
	m_UsePaddingValue = false;
//	m_UseAllPixels = false;
	m_DerivativeStepLength = 0.1;
	m_DerivativeStepLengthScales.Fill(1);
	m_UpperBoundIncreaseFactor = 0.001;
	m_PaddingValue = NumericTraits<FixedImagePixelType>::Zero;
	m_Histogram = HistogramType::New();
#ifdef ITK_USE_REVIEW_STATISTICS
	m_Histogram->SetMeasurementVectorSize(2);
	m_Thresholds.SetSize(2);
	m_BinThresholds.SetSize(2);
#endif
	m_LowerBoundSetByUser = false;
	m_UpperBoundSetByUser = false;
	m_Max.resize(2);
	m_Min.resize(2);

	m_ThresholdsSetByUser = false;
	m_WeightingFactor = 0.0;
	m_UseCoincidenceWeighting = false;
	m_NumberOfSpatialSamples = 0;
	m_BinThresholds[0] = 0;
	m_BinThresholds[1] = 0;
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::SetUpperBound(
		const MeasurementVectorType & bounds) {
	m_UpperBound = bounds;
	m_UpperBoundSetByUser = true;
	this->Modified();
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::SetLowerBound(
		const MeasurementVectorType & bounds) {
	m_LowerBound = bounds;
	m_LowerBoundSetByUser = true;
	this->Modified();
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::SetBackgroundThresholds(
		const MeasurementVectorType & thres) {
	m_Thresholds = thres;
	m_ThresholdsSetByUser = true;
	this->Modified();
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::Initialize()
		throw (ExceptionObject) {
	Superclass::Initialize();

	if (!this->m_FixedImage) {
		itkExceptionMacro(<<"Fixed image has not been set.");
	} else if (!this->m_MovingImage) {
		itkExceptionMacro(<<"Moving image has not been set.");
	}

	if (!m_LowerBoundSetByUser || !m_UpperBoundSetByUser) {
		// Calculate min and max image values in fixed image.
		FixedImageConstPointerType pFixedImage = this->m_FixedImage;
		ImageRegionConstIterator<FixedImageType> fiIt(pFixedImage,
				pFixedImage-> GetBufferedRegion());
		fiIt.GoToBegin();
		FixedImagePixelType minFixed = fiIt.Value();
		FixedImagePixelType maxFixed = fiIt.Value();
		++fiIt;
		while (!fiIt.IsAtEnd()) {
			FixedImagePixelType value = fiIt.Value();

			if (value < minFixed) {
				minFixed = value;
			} else if (value > maxFixed) {
				maxFixed = value;
			}

			++fiIt;
		}

		// Calculate min and max image values in moving image.
		MovingImageConstPointerType pMovingImage = this->m_MovingImage;
		ImageRegionConstIterator<MovingImageType> miIt(pMovingImage,
				pMovingImage-> GetBufferedRegion());
		miIt.GoToBegin();
		MovingImagePixelType minMoving = miIt.Value();
		MovingImagePixelType maxMoving = miIt.Value();
		++miIt;
		while (!miIt.IsAtEnd()) {
			MovingImagePixelType value = miIt.Value();

			if (value < minMoving) {
				minMoving = value;
			} else if (value > maxMoving) {
				maxMoving = value;
			}
			++miIt;
		}

		// Initialize the upper and lower bounds of the histogram.
		if (!m_LowerBoundSetByUser) {
#ifdef ITK_USE_REVIEW_STATISTICS
			m_LowerBound.SetSize(2);
#endif
			m_LowerBound[0] = minFixed;
			m_LowerBound[1] = minMoving;
		}

		if (!m_UpperBoundSetByUser) {
#ifdef ITK_USE_REVIEW_STATISTICS
			m_UpperBound.SetSize(2);
#endif
			m_UpperBound[0] = maxFixed + (maxFixed - minFixed)
					* m_UpperBoundIncreaseFactor;
			m_UpperBound[1] = maxMoving + (maxMoving - minMoving)
					* m_UpperBoundIncreaseFactor;
		}

	}

	HistogramSizeType fixedSize;
	fixedSize.SetSize(1);
	fixedSize.Fill(m_HistogramSize[0]);
	HistogramSizeType movingSize;
	movingSize.SetSize(1);
	movingSize.Fill(m_HistogramSize[1]);

	if (m_UseAdaptativeBining) {
		// TODO: Should I use masks???
		typename FixedAdaptorType::Pointer fixedAdaptor =
				FixedAdaptorType::New();
		fixedAdaptor->SetImage(this->m_FixedImage);
		typename FixedHistogramFilter::Pointer fixedHistFilter =
				FixedHistogramFilter::New();
		fixedHistFilter->SetInput(fixedAdaptor);
		fixedHistFilter->SetHistogramSize(fixedSize);
		fixedHistFilter->SetMaxEstimatorIterations(7);
		fixedHistFilter->SetMarginalScale(10);
		//  fixedHistFilter->SetHistogramBinMinimum( m_LowerBound[0] );
		//  fixedHistFilter->SetHistogramBinMaximum( m_UpperBound[0] );
		//  fixedHistFilter->SetAutoMinimumMaximum( false );
		fixedHistFilter->Update();

		const HistogramType * fixedHistogram = fixedHistFilter->GetOutput();
		typename MovingAdaptorType::Pointer movingAdaptor =
				MovingAdaptorType::New();
		movingAdaptor->SetImage(this->m_MovingImage);
		typename MovingHistogramFilter::Pointer movingHistFilter =
				MovingHistogramFilter::New();
		movingHistFilter->SetInput(movingAdaptor);
		movingHistFilter->SetHistogramSize(movingSize);
		movingHistFilter->SetMaxEstimatorIterations(7);
		movingHistFilter->SetMarginalScale(10);
		//movingHistFilter->SetHistogramBinMinimum( m_LowerBound[1] );
		//movingHistFilter->SetHistogramBinMaximum( m_UpperBound[1] );
		//movingHistFilter->SetAutoMinimumMaximum( false );
		movingHistFilter->Update();

		const HistogramType * movingHistogram = movingHistFilter->GetOutput();

		m_Max[0] = fixedHistogram->GetDimensionMaxs(0);
		m_Min[0] = fixedHistogram->GetDimensionMins(0);
		m_Max[1] = movingHistogram->GetDimensionMaxs(0);
		m_Min[1] = movingHistogram->GetDimensionMins(0);

		if (m_UseCoincidenceWeighting && !m_ThresholdsSetByUser) {
			fixedSize.Fill(2);
			movingSize.Fill(2);

			fixedHistFilter->SetHistogramSize(fixedSize);
			movingHistFilter->SetHistogramSize(movingSize);
			fixedHistFilter->Update();
			movingHistFilter->Update();

			m_Thresholds[0] = fixedHistFilter->GetOutput()->GetBinMax(0, 0);
			m_Thresholds[1] = movingHistFilter->GetOutput()->GetBinMax(0, 0);

			itkDebugMacro("Coincidence Weighting->Thresholds Computed = " << this->m_Thresholds );
		}
	}
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::SetTransform(
		TransformType * transform) {
	if (m_DerivativeStepLengthScales.GetSize()
			!= transform->GetNumberOfParameters()) {
		m_DerivativeStepLengthScales.SetSize(transform->GetNumberOfParameters());
		m_DerivativeStepLengthScales.Fill(1.0);
	}
	Superclass::SetTransform(transform);
}

template<class TFixedImage, class TMovingImage>
typename HistogramImageToImageMetric<TFixedImage, TMovingImage>::MeasureType HistogramImageToImageMetric<
		TFixedImage, TMovingImage>::GetValue(
		const TransformParametersType& parameters) const {
	itkDebugMacro("GetValue( " << parameters << " ) ");

	this->ComputeHistogram(parameters, *m_Histogram);

	if (m_UseAdaptativeBining && m_UseCoincidenceWeighting && m_WeightingFactor
			== 0.0) {
		double pB = 0.0;
		double pO = 0.0;

		typename HistogramType::MeasurementVectorType actualSample;
#ifdef ITK_USE_REVIEW_STATISTICS
		actualSample.SetSize(2);
#endif

		typename HistogramType::ConstIterator it = m_Histogram->Begin();
		typename HistogramType::ConstIterator end = m_Histogram->End();

		while (it != end) {
			unsigned int freq = it.GetFrequency();
			actualSample = it.GetMeasurementVector();

			if (actualSample[0] <= m_Thresholds[0] && actualSample[1]
					<= m_Thresholds[1]) {
				// Add to pB
				pB += freq;
			} else {
				// Add to pO
				pO += freq;
			}
			++it;
		}

		pB /= m_Histogram->GetTotalFrequency();
		pO /= m_Histogram->GetTotalFrequency();

		if (pB + pO < 0.98 || pB + pO > 1.02) {
			itkExceptionMacro(
					<< "Object and Background Probabilities are not complementary");
		}

		double ratio = pO / pB;
		m_WeightingFactor = (ratio < 1.0) ? ratio : 1.0;
		itkDebugMacro("Weighting Factor = " << m_WeightingFactor );

		this->ComputeHistogram(parameters, *m_Histogram);
	}

	return this->EvaluateMeasure(*m_Histogram);
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::GetDerivative(
		const TransformParametersType& parameters, DerivativeType& derivative) const {
	itkDebugMacro("GetDerivative( " << parameters << " ) ");

	const unsigned int ParametersDimension = this->GetNumberOfParameters();

	// Make sure the scales have been set
	if (m_DerivativeStepLengthScales.size() != ParametersDimension) {
		itkExceptionMacro(<< "The size of DerivativesStepLengthScales is "
				<< m_DerivativeStepLengthScales.size()
				<< ", but the Number of Parameters is "
				<< ParametersDimension
				<< ".");
	}

	// Calculate gradient.
	derivative = DerivativeType(ParametersDimension);
	derivative.Fill(NumericTraits<ITK_TYPENAME
	DerivativeType::ValueType>::Zero);

	typename HistogramType::Pointer pHistogram = HistogramType::New();
#ifdef ITK_USE_REVIEW_STATISTICS
	pHistogram->SetMeasurementVectorSize(2);
#endif
	this->ComputeHistogram(parameters, *pHistogram);

	TransformParametersType newParameters1;
	TransformParametersType newParameters2;

	for (unsigned int i = 0; i < ParametersDimension; i++) {
		typename HistogramType::Pointer pHistogram2 = HistogramType::New();
#ifdef ITK_USE_REVIEW_STATISTICS
		pHistogram2->SetMeasurementVectorSize(2);
#endif
		this->CopyHistogram(*pHistogram2, *pHistogram);

		TransformParametersType newParameters = parameters;
		newParameters[i] -= m_DerivativeStepLength
				/ m_DerivativeStepLengthScales[i];
		this->ComputeHistogram(newParameters, *pHistogram2);

		MeasureType e0 = EvaluateMeasure(*pHistogram2);

		pHistogram2 = HistogramType::New();
#ifdef ITK_USE_REVIEW_STATISTICS
		pHistogram2->SetMeasurementVectorSize(2);
#endif
		this->CopyHistogram(*pHistogram2, *pHistogram);

		newParameters = parameters;
		newParameters[i] += m_DerivativeStepLength
				/ m_DerivativeStepLengthScales[i];

		this->ComputeHistogram(newParameters, *pHistogram2);

		MeasureType e1 = EvaluateMeasure(*pHistogram2);

		derivative[i] = (e1 - e0) / (2 * m_DerivativeStepLength
				/ m_DerivativeStepLengthScales[i]);
	}

	//std::cout << derivative << std::endl;
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::GetValueAndDerivative(
		const TransformParametersType& parameters, MeasureType& value,
		DerivativeType& derivative) const {
	value = GetValue(parameters);
	this->GetDerivative(parameters, derivative);
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::ComputeHistogram(
		TransformParametersType const& parameters, HistogramType& histogram) const {
	FixedImageConstPointerType fixedImage = this->m_FixedImage;

	if (!fixedImage) {
		itkExceptionMacro(<< "Fixed image has not been assigned");
	}
	this->m_NumberOfPixelsCounted = 0;
	this->SetTransformParameters(parameters);

	typedef std::vector<typename HistogramType::BinMaxVectorType>
			BinMaxContainerType;
	typedef std::vector<typename HistogramType::BinMinVectorType>
			BinMinContainerType;
	if (m_UseAdaptativeBining) {
		histogram.Initialize((BinMaxContainerType) m_Min,
				(BinMinContainerType) m_Max);
		m_BinThresholds = m_Histogram->GetIndex(m_Thresholds);
	} else {
		histogram.Initialize(m_HistogramSize, m_LowerBound, m_UpperBound);
	}


	if (this->m_UseSequentialSampling ) {
		typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType>
				FixedIteratorType;

		typename FixedImageType::IndexType index;
		typename FixedImageType::RegionType fixedRegion;

		fixedRegion = this->GetFixedImageRegion();
		FixedIteratorType ti(fixedImage, fixedRegion);

		ti.GoToBegin();
		while (!ti.IsAtEnd()) {
			index = ti.GetIndex();

			if (fixedRegion.IsInside(index) && (!m_UsePaddingValue
					|| (m_UsePaddingValue && ti.Get() > m_PaddingValue))) {
				InputPointType inputPoint;
				fixedImage->TransformIndexToPhysicalPoint(index, inputPoint);

				if (this->m_FixedImageMask
						&& !this->m_FixedImageMask->IsInside(inputPoint)) {
					++ti;
					continue;
				}

				OutputPointType transformedPoint =
						this->m_Transform->TransformPoint(inputPoint);

				if (this->m_MovingImageMask
						&& !this->m_MovingImageMask->IsInside(transformedPoint)) {
					++ti;
					continue;
				}

				if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
					const RealType movingValue =
							this->m_Interpolator->Evaluate(transformedPoint);
					const RealType fixedValue = ti.Get();
					this->m_NumberOfPixelsCounted++;

					typename HistogramType::MeasurementVectorType sample;
#ifdef ITK_USE_REVIEW_STATISTICS
					sample.SetSize(2);
#endif
					sample[0] = fixedValue;
					sample[1] = movingValue;

					unsigned int freq = 1;

					if (m_UseCoincidenceWeighting && m_WeightingFactor > 0.0
							&& m_WeightingFactor < 1.0) {
						freq = 100;
						if (fixedValue < m_Thresholds[0] && movingValue
								< m_Thresholds[1]) {
							double D = sqrt(pow((double) m_Thresholds[0],
									(int) 2) + pow((double) m_Thresholds[1],
									(int) 2));
							double distance = sqrt(pow((double) fixedValue,
									(int) 2) + pow((double) movingValue,
									(int) 2)) / D;

							if (distance > 1.0)
								distance = 1.0;

							freq *= distance * m_WeightingFactor;

							if (freq > (m_WeightingFactor * 100)) {
								std::cout << "error" << std::endl;
							}
						}
					}
					histogram.IncreaseFrequency(sample, freq);
				}
			}

			++ti;
		}
	} else {
		this->FillSubsampledHistogram(histogram);
	}

	itkDebugMacro("NumberOfPixelsCounted = " << this->m_NumberOfPixelsCounted );
	if (this->m_NumberOfPixelsCounted == 0) {
		itkExceptionMacro(
				<< "All the points mapped to outside of the moving image");
	}
	/*
	 if ( m_UseCoincidenceWeighting && m_WeightingFactor < 1.0 && m_WeightingFactor > 0.0 )
	 {
	 double D = sqrt(pow((double) m_BinThresholds[0], (int) 2) + pow( (double) m_BinThresholds[1], (int) 2));
	 for ( unsigned int i = 0; i < m_BinThresholds[0]; i++ )
	 {
	 for ( unsigned int j = 0; j < m_BinThresholds[1]; j++ )
	 {
	 typename HistogramType::IndexType index;
	 index.SetSize(2);
	 index[0] = i;
	 index[1] = j;
	 double distance = sqrt(pow( (double) i,(int) 2) + pow( (double) j,(int) 2)) / D;
	 unsigned int freq = histogram.GetFrequency(index);
	 histogram.SetFrequency( index, (unsigned int) freq*m_WeightingFactor * distance );
	 }
	 }
	 }
	 */
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::FillSubsampledHistogram(
		HistogramType& histogram) const {
	typedef itk::ImageRandomConstIteratorWithIndex<FixedImageType>
			FixedIteratorType;

	typename FixedImageType::IndexType index;
	typename FixedImageType::RegionType fixedRegion;

	fixedRegion = this->GetFixedImageRegion();
	FixedIteratorType ti(this->m_FixedImage, fixedRegion);

	this->m_NumberOfPixelsCounted = 0;

	ti.SetNumberOfSamples(m_NumberOfSpatialSamples);
	ti.GoToBegin();
	while (!ti.IsAtEnd()) {
		index = ti.GetIndex();

		if (fixedRegion.IsInside(index) && (!m_UsePaddingValue
				|| (m_UsePaddingValue && ti.Get() > m_PaddingValue))) {
			InputPointType inputPoint;
			this->m_FixedImage->TransformIndexToPhysicalPoint(index, inputPoint);

			if (this->m_FixedImageMask && !this->m_FixedImageMask->IsInside(
					inputPoint)) {
				++ti;
				continue;
			}

			OutputPointType transformedPoint =
					this->m_Transform->TransformPoint(inputPoint);

			if (this->m_MovingImageMask && !this->m_MovingImageMask->IsInside(
					transformedPoint)) {
				++ti;
				continue;
			}

			if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
				const RealType movingValue = this->m_Interpolator->Evaluate(
						transformedPoint);
				const RealType fixedValue = ti.Get();
				this->m_NumberOfPixelsCounted++;

				typename HistogramType::MeasurementVectorType sample;
#ifdef ITK_USE_REVIEW_STATISTICS
				sample.SetSize(2);
#endif
				sample[0] = fixedValue;
				sample[1] = movingValue;

				unsigned int freq = 1.0;

				if (m_UseCoincidenceWeighting && m_WeightingFactor > 0.0
						&& m_WeightingFactor < 1.0) {
					freq = 100;
					HistogramType::IndexType idx = histogram.GetIndex(sample);

					if (idx[0] <= m_BinThresholds[0] && idx[1]
							<= m_BinThresholds[1]) {
						//double D = sqrt(pow( (double) m_BinThresholds[0],(int) 2) + pow( (double) m_BinThresholds[1],(int) 2));
						//double distance = sqrt(pow( (double) idx[0],(int) 2) + pow( (double) idx[1],(int) 2)) / D;
						freq *= m_WeightingFactor;
					}
				}

				histogram.IncreaseFrequency(sample, freq);
			}
		}

		++ti;
	}

	//   HistogramType::IndexType idx;
	//#ifdef ITK_USE_REVIEW_STATISTICS
	//	idx.SetSize(2);
	//#endif
	//idx[0] = 0;
	//idx[1] = 0;


	//histogram.SetFrequency( idx, 0 );
	//idx[0]++;
	//histogram.SetFrequency( idx, 0 );
	//idx[0] = 0;
	//idx[1] = 1;
	//histogram.SetFrequency( idx, 0 );

	if (m_UseCoincidenceWeighting && m_WeightingFactor > 0.0
			&& m_WeightingFactor < 1.0) {
		HistogramType::Iterator it = histogram.Begin();
		HistogramType::Iterator end = histogram.End();

		it.SetFrequency(0);
		++it;

		while (it != end) {
			it.SetFrequency(it.GetFrequency() / 100.0);
			++it;
		}

	}
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::CopyHistogram(
		HistogramType& target, HistogramType& source) const {
	// Initialize the target.


	typename HistogramType::SizeType size = source.GetSize();

	if (m_UseAdaptativeBining) {
		typedef std::vector<typename HistogramType::BinMaxVectorType>
				BinMaxContainerType;
		typedef std::vector<typename HistogramType::BinMinVectorType>
				BinMinContainerType;
		BinMaxContainerType maxs;
		BinMinContainerType mins;

		for (unsigned int dim = 0; dim < size.Size(); dim++) {
			maxs.push_back(source.GetDimensionMins(dim));
			mins.push_back(source.GetDimensionMins(dim));
		}

		target.Initialize(mins, maxs);
	} else {
		typename HistogramType::MeasurementVectorType min, max;
#ifdef ITK_USE_REVIEW_STATISTICS
		min.SetSize(2);
		max.SetSize(2);
#endif
		for (unsigned int i = 0; i < min.Size(); i++) {
			min[i] = source.GetBinMin(i, 0);
		}

		for (unsigned int i = 0; i < max.Size(); i++) {
			max[i] = source.GetBinMax(i, size[i] - 1);
		}

		target.Initialize(size, min, max);
	}

	// Copy the values.
	typename HistogramType::Iterator sourceIt = source.Begin();
	typename HistogramType::Iterator sourceEnd = source.End();
	typename HistogramType::Iterator targetIt = target.Begin();
	typename HistogramType::Iterator targetEnd = target.End();

	while (sourceIt != sourceEnd && targetIt != targetEnd) {
#ifdef ITK_USE_REVIEW_STATISTICS
		typename HistogramType::AbsoluteFrequencyType
#else
		typename HistogramType::FrequencyType
#endif
		freq = sourceIt.GetFrequency();

		if (freq > 0) {
			targetIt.SetFrequency(freq);
		}

		++sourceIt;
		++targetIt;
	}
}

template<class TFixedImage, class TMovingImage>
void HistogramImageToImageMetric<TFixedImage, TMovingImage>::PrintSelf(
		std::ostream& os, Indent indent) const {
	Superclass::PrintSelf(os, indent);
	os << indent << "Use Coincidence Weighting?: " << m_UseCoincidenceWeighting
			<< std::endl;
	if (m_UseCoincidenceWeighting) {
		os << indent << indent << "* Weighting Factor: " << m_WeightingFactor
				<< std::endl;
		os << indent << indent << "* Thresholds (#bin): ["
				<< m_BinThresholds[0] << ", " << m_BinThresholds[1] << "]"
				<< std::endl;
		os << indent << indent << "* Thresholds (values): [" << m_Thresholds[0]
				<< ", " << m_Thresholds[1] << "]" << std::endl;
	}
	os << indent << "Use Adaptative Bining?: " << m_UseAdaptativeBining
			<< std::endl;
	os << indent << "Padding value: " << static_cast<typename NumericTraits<
			FixedImagePixelType>::PrintType> (m_PaddingValue) << std::endl;
	os << indent << "Use padding value?: " << m_UsePaddingValue << std::endl;
	os << indent << "Derivative step length: " << m_DerivativeStepLength
			<< std::endl;
	os << indent << "Derivative step length scales: ";
	os << m_DerivativeStepLengthScales << std::endl;
	os << indent << "Histogram size: ";
	os << m_HistogramSize << std::endl;
	os << indent << "Histogram upper bound increase factor: ";
	os << m_UpperBoundIncreaseFactor << std::endl;
	os << indent << "Histogram computed by GetValue(): ";
	os << m_Histogram.GetPointer() << std::endl;
}

} // end namespace itk

#endif // itkHistogramImageToImageMetric_txx
