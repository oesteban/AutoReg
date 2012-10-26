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
#ifndef METRICHELPER_TXX
#define METRICHELPER_TXX

#ifndef DEG
#define DEG M_PI/180.0
#endif

#ifndef NUMBER_OF_THREADS
#define NUMBER_OF_THREADS 4
#endif

#include "MetricHelper.h"
template<class TFixedImage, class TMovingImage>
MetricHelper<TFixedImage, TMovingImage>::MetricHelper() {
	SetMetricType(MI);
}

template<class TFixedImage, class TMovingImage>
void MetricHelper<TFixedImage, TMovingImage>::SetMetricType(METRIC_TYPE type) {
	m_Type = type;

	// Set default config
	typename MIType::HistogramType::SizeType histogramSize;
	unsigned long histogramBins = 70;
	histogramSize.SetSize(2);
	histogramSize.Fill(histogramBins);

	this->SetDerivativeStepLength(1.0);

	m_Scales = ScalesType(6);
	m_Scales[0] = m_Scales[1] = m_Scales[2] = 1.0 / (DEG * 0.1); // DEG*4
	m_Scales[3] = m_Scales[4] = m_Scales[5] = 1.0; // 1.0 / 4.0

	this->SetDerivativeStepLengthScales(m_Scales);

	switch (m_Type) {
	case MI:
	case NMI:
		m_Parameters.resize(9);
		m_Parameters[0] = histogramSize;
		m_Parameters[4] = false;
		m_Parameters[5] = false;
		m_Parameters[6] = 50000u;
		m_Parameters[7] = NUMBER_OF_THREADS;
		m_Parameters[8] = false;
		break;
	case MATTES:
		m_Parameters.resize(6);
		m_Parameters[0] = histogramBins;
		m_Parameters[1] = 20000u;
		m_Parameters[2] = true;
		m_Parameters[3] = 123u;
		m_Parameters[4] = NUMBER_OF_THREADS;
		m_Parameters[5] = -255.0;
		break;
	case VIOLA_WELLS:
		m_Parameters.resize(0);
		break;
	default:
		break;
	}

}

template<class TFixedImage, class TMovingImage>
void MetricHelper<TFixedImage, TMovingImage>::SetFixedImageMask( const typename MetricHelper<TFixedImage, TMovingImage>::FixedImageType* fixedMask ) {
	typename FixedMaskCasterType::Pointer cast = FixedMaskCasterType::New();
	cast->SetInput( fixedMask );
	cast->Update();
	m_FixedMask = FixedMaskSpatialObjectType::New();
	typename FixedMaskCasterType::OutputImageType::Pointer im = cast->GetOutput();
	m_FixedMask->SetImage( im );
}

template<class TFixedImage, class TMovingImage>
void MetricHelper<TFixedImage, TMovingImage>::SetMovingImageMask( const typename MetricHelper<TFixedImage, TMovingImage>::MovingImageType* movingMask ){
	typename MovingMaskCasterType::Pointer cast = MovingMaskCasterType::New();
	cast->SetInput( movingMask );
	cast->Update();
	m_MovingMask = MovingMaskSpatialObjectType::New();
	m_MovingMask->SetImage( cast->GetOutput() );
}

template<class TFixedImage, class TMovingImage>
template<class TRegistration>
void MetricHelper<TFixedImage, TMovingImage>::Connect(TRegistration* r) {
	switch (m_Type) {
	case NMI:
		m_NMIMetric = NMIType::New();
		try {
			m_NMIMetric->SetHistogramSize(boost::any_cast<
					typename NMIType::HistogramType::SizeType>(m_Parameters[0]));
		} catch (...) {
		}

		try {
			m_NMIMetric->SetLowerBound(boost::any_cast<
					typename NMIType::MeasurementVectorType>(m_Parameters[1]));
		} catch (...) {
		}

		try {
			m_NMIMetric->SetUpperBound(boost::any_cast<
					typename NMIType::MeasurementVectorType>(m_Parameters[2]));
		} catch (...) {
		}

		try {
			m_NMIMetric->SetPaddingValue(boost::any_cast<double>(
					m_Parameters[3]));
			m_NMIMetric->SetUsePaddingValue(true);
		} catch (...) {
		}
#ifdef USE_ADAPTATIVE_BINING
		try {
			m_NMIMetric->SetUseCoincidenceWeighting( boost::any_cast< bool > ( m_Parameters[4] ) );
		} catch (...) {}
		try {
			m_NMIMetric->SetUseAdaptativeBining( boost::any_cast< bool > ( m_Parameters[5] ) );
		} catch (...) {}
		try {
			m_NMIMetric->SetNumberOfSpatialSamples( boost::any_cast< unsigned int > ( m_Parameters[6] ) );
		} catch (...) {}
		try {
			m_NMIMetric->SetUseAllPixels( boost::any_cast< bool > ( m_Parameters[8] ) );
		} catch (...) {}
#endif
		try {
			m_NMIMetric->SetNumberOfThreads(boost::any_cast<unsigned int >(
					m_Parameters[7]));
		} catch (...) {
		}

		if (m_StepLength > 0)
			m_NMIMetric->SetDerivativeStepLength(m_StepLength);
		if (!m_Scales.empty())
			m_NMIMetric->SetDerivativeStepLengthScales(m_Scales);

		m_NMIMetric->SetFixedImageMask( m_FixedMask );
		m_NMIMetric->SetMovingImageMask( m_MovingMask );

		r->SetMetric(m_NMIMetric);
		break;
	case MATTES:
		m_MattesMetric = MattesType::New();
		try {
			m_MattesMetric->SetNumberOfHistogramBins(boost::any_cast<
					unsigned long >(m_Parameters[0]));
		} catch (...) {
		}

		try {
			m_MattesMetric->SetNumberOfSpatialSamples(boost::any_cast<
					unsigned int >(m_Parameters[1]));
		} catch (...) {
		}

		try {
			m_MattesMetric->SetUseExplicitPDFDerivatives(boost::any_cast<bool>(
					m_Parameters[2]));
		} catch (...) {
		}

		try {
			m_MattesMetric->ReinitializeSeed(boost::any_cast<unsigned int >(
					m_Parameters[3]));
		} catch (...) {
		}

		try {
			m_MattesMetric->SetNumberOfThreads(boost::any_cast<unsigned int >(
					m_Parameters[4]));
		} catch (...) {
		}

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
		try {
			double thres = boost::any_cast<double>(m_Parameters[5]);
			if (thres > -255.0) {
				m_MattesMetric->SetFixedImageSamplesIntensityThreshold(thres);
				m_MattesMetric->SetUseFixedImageSamplesIntensityThreshold(true);
			}
		} catch (std::exception& e) {
		}
#endif

		m_MattesMetric->SetFixedImageMask( m_FixedMask );
		m_MattesMetric->SetMovingImageMask( m_MovingMask );

		r->SetMetric(m_MattesMetric);
		break;
	case VIOLA_WELLS:
		m_ViolaWellsMetric = ViolaWellsType::New();
		try {
			//	m_ViolaWellsMetric->SetHistogramSize( boost::any_cast< ViolaWellsType::HistogramType::SizeType> ( m_Parameters[0] ) );
		} catch (...) {
		}
		m_ViolaWellsMetric->SetFixedImageMask( m_FixedMask );
		m_ViolaWellsMetric->SetMovingImageMask( m_MovingMask );
		r->SetMetric(m_ViolaWellsMetric);
		break;

	case MI: // TODO update parameters!
	default:
		m_MIMetric = MIType::New();
		try {
			m_MIMetric->SetHistogramSize(boost::any_cast<
					typename MIType::HistogramType::SizeType>(m_Parameters[0]));
		} catch (...) {
		}
#ifdef USE_ADAPTATIVE_BINING
		try {
			m_MIMetric->SetUseCoincidenceWeighting( boost::any_cast< bool > ( m_Parameters[1] ) );
		} catch (...) {}
		try {
			m_MIMetric->SetUseAdaptativeBining( boost::any_cast< bool > ( m_Parameters[2] ) );
		} catch (...) {}
		try {
			m_MIMetric->SetNumberOfSpatialSamples( boost::any_cast< unsigned int > ( m_Parameters[3] ) );
		} catch (...) {}
#endif
		try {
			m_MIMetric->SetNumberOfThreads(boost::any_cast<unsigned int >(
					m_Parameters[4]));
		} catch (...) {
		}
		m_MIMetric->SetFixedImageMask( m_FixedMask );
		m_MIMetric->SetMovingImageMask( m_MovingMask );
		r->SetMetric(m_MIMetric);
		break;
	}

}
;

template<class TFixedImage, class TMovingImage>
unsigned int MetricHelper<TFixedImage, TMovingImage>::GetNumberOfParameters() {
	return m_Parameters.size();
}
;

template<class TFixedImage, class TMovingImage>
template<class TMetric>
TMetric* MetricHelper<TFixedImage, TMovingImage>::GetMetric() {
	switch (m_Type) {
	case MI:
		return m_MIMetric;
	case NMI:
		return m_NMIMetric;
	case MATTES:
		return m_MattesMetric;
	case VIOLA_WELLS:
		return m_ViolaWellsMetric;
	default:
		return NULL;
	}
}

template<class TFixedImage, class TMovingImage>
void MetricHelper<TFixedImage, TMovingImage>::PrintSelf(std::ostream& os,
		itk::Indent indent) const {

	Superclass::PrintSelf(os, indent);

	switch (m_Type) {
	case MI:
		m_MIMetric->Print(os, indent);
		break;
	case NMI:
		if (m_NMIMetric.IsNotNull())
			m_NMIMetric->Print(os, indent);
		break;
	case MATTES:
		m_MattesMetric->Print(os, indent);
		break;
	case VIOLA_WELLS:
		m_ViolaWellsMetric->Print(os, indent);
		break;
	default:
		os << indent << "No metric configured" << std::endl;
	}

}
;

#endif // METRICHELPER_TXX
