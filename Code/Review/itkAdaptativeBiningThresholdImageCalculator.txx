/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkAdaptativeBiningThresholdImageCalculator.txx,v $
 Language:  C++
 Date:      $Date: 2009-01-26 21:45:54 $
 Version:   $Revision: 1.9 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkAdaptativeBiningThresholdImageCalculator_txx
#define __itkAdaptativeBiningThresholdImageCalculator_txx

#include "itkAdaptativeBiningThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk {

/**
 * Constructor
 */
template<class TInputImage>
AdaptativeBiningThresholdImageCalculator<TInputImage>::AdaptativeBiningThresholdImageCalculator() {
	itkDebugMacro("Constructor");
	m_Image = NULL;
	m_Threshold = NumericTraits<PixelType>::Zero;
	m_NumberOfPDFBins = 128;
	m_RegionSetByUser = false;

	m_PDF = HistogramType::New();

	#ifdef ITK_USE_REVIEW_STATISTICS
	  m_PDF->SetMeasurementVectorSize(1);
	#endif


}

/*
 * Compute the Otsu's threshold
 */
template<class TInputImage>
void AdaptativeBiningThresholdImageCalculator<TInputImage>::Compute(void) {
	if (!m_Image) {
		return;
	}

	// TODO manage regions
	if (!m_RegionSetByUser) {
		m_Region = m_Image->GetRequestedRegion();
	}

	// Generate histogram
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage( m_Image );

	typename HistogramFilter::Pointer histFilter = HistogramFilter::New();
	histFilter->SetInput( adaptor );

	typename HistogramFilter::HistogramSizeType s;
	s.SetSize( 1 );
	s.Fill( 2 );
	histFilter->SetHistogramSize( s );
	histFilter->SetMaxEstimatorIterations( 7 );
	histFilter->SetMarginalScale( 10 );
	histFilter->SetPDFDensityFactor( 56 );
	histFilter->Update();

	const HistogramType * histogram = histFilter->GetOutput();

	m_Threshold = histogram->GetBinMax(0,0);
}

template<class TInputImage>
void AdaptativeBiningThresholdImageCalculator<TInputImage>::SetRegion(
		const RegionType & region) {
	m_Region = region;
	m_RegionSetByUser = true;
}

template<class TInputImage>
void AdaptativeBiningThresholdImageCalculator<TInputImage>::PrintSelf(
		std::ostream& os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Threshold: " << m_Threshold << std::endl;
	os << indent << "NumberOfHistogramBins: " << m_NumberOfPDFBins << std::endl;
	os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
