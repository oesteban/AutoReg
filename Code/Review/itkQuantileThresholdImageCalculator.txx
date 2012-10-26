/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkQuantileThresholdImageCalculator.txx,v $
 Language:  C++
 Date:      $Date: 2009-01-26 21:45:54 $
 Version:   $Revision: 1.9 $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkQuantileThresholdImageCalculator_txx
#define __itkQuantileThresholdImageCalculator_txx

#include "itkQuantileThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk {

/**
 * Constructor
 */
template<class TInputImage>
QuantileThresholdImageCalculator<TInputImage>::QuantileThresholdImageCalculator()
{
	itkDebugMacro("Constructor");
    m_Percentage = 0.0;
	m_Image=NULL;
	m_NumberOfHistogramBins=128;
	m_Threshold=0.0;
	m_HistogramComputed=false;
}

/*
 * Compute the Otsu's threshold
 */
template<class TInputImage>
void QuantileThresholdImageCalculator<TInputImage>::Compute(void) {
	if (!m_Image) {
		return;
    }

  if (!m_HistogramComputed)
  {
	// Generate histogram
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	adaptor->SetImage( m_Image );

	typename HistogramFilter::Pointer histFilter = HistogramFilter::New();
	histFilter->SetInput( adaptor );

	typename HistogramFilter::HistogramSizeType s;
	s.SetSize( 1 );
	s.Fill( m_NumberOfHistogramBins );
	histFilter->SetHistogramSize( s );
	histFilter->Update();

	m_Histogram = histFilter->GetOutput();
    m_HistogramComputed=true;
  }
	m_Threshold = m_Histogram->Quantile(0, m_Percentage);
}

template<class TInputImage>
void QuantileThresholdImageCalculator<TInputImage>::PrintSelf(
		std::ostream& os, Indent indent) const {
	Superclass::PrintSelf(os, indent);

	os << indent << "Threshold: " << m_Threshold << std::endl;
	os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
	os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
