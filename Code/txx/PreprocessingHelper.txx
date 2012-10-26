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

    File: PreprocessingHelper.txx
	Project: AutoRegistration
    Created on: 18/11/2010
    Author: oesteban
 
*/

#ifndef PREPROCESSINGHELPER_TXX_
#define PREPROCESSINGHELPER_TXX_

#include <PreprocessingHelper.h>


template <class TInputImage, class TOutputImage>
const float PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_POW   =  2.0;

template <class TInputImage, class TOutputImage>
const float PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_SIGMA = 10.0;

template <class TInputImage, class TOutputImage>
PreprocessingHelper<TInputImage,TOutputImage>::PreprocessingHelper(): m_QuantileComputed(false), m_RATSComputed(false) {
   m_MinThreshold = itk::NumericTraits<InputPixelType>::min();
   m_Min = itk::NumericTraits<InputPixelType>::min();

   m_MaxThreshold = itk::NumericTraits<InputPixelType>::max();
   m_Max = itk::NumericTraits<InputPixelType>::max();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::ComputeQuantileThresholds(float minPercentage, float maxPercentage )
{
	typename QuantileCalculator::Pointer th_calc = QuantileCalculator::New();
	th_calc->SetImage(m_InputImage);
	th_calc->SetPercentage( minPercentage );
	th_calc->SetNumberOfHistogramBins(256);
	th_calc->Compute();

	m_MinThreshold = th_calc->GetThreshold();

	th_calc->SetPercentage( maxPercentage );
	th_calc->Compute();

	m_MaxThreshold = th_calc->GetThreshold();

	m_QuantileComputed=true;
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::ComputeRATSThreshold(float pow, float sigma)
{
	typename GradientType::Pointer gradient = GradientType::New();
	gradient->SetInput( m_InputImage );
	gradient->SetSigma( sigma );
	gradient->Update();

	// Compute the Threshold for the input image
	typename CalculatorType::Pointer thresholdCalculator = CalculatorType::New();
	thresholdCalculator->SetInput( m_InputImage );
	thresholdCalculator->SetGradient( gradient->GetOutput() );
	thresholdCalculator->SetPow( pow );
	thresholdCalculator->Compute();

	m_RATSThreshold = thresholdCalculator->GetOutput();
	m_RATSComputed = true;
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetHistogramMatchingImage(const typename PreprocessingHelper<TInputImage,TOutputImage>::OutputImageType * refImage) {
	typename HistogramFilter::Pointer f = HistogramFilter::New();
	f->SetInput( m_OutputImage );
	f->SetReferenceImage( refImage );
	f->SetNumberOfHistogramLevels( 128 );
	f->SetNumberOfMatchPoints( 56 );
	f->ThresholdAtMeanIntensityOn();
	f->Update();
	m_OutputImage = f->GetOutput();

	typename OutCalculator::Pointer stats = OutCalculator::New();
	stats->SetInput( m_OutputImage );
	stats->Update();

	m_Max = stats->GetMaximum();
	m_Min = stats->GetMinimum();
	m_Mean= stats->GetMean();
	m_Sigma = stats->GetSigma();

}

template <class TInputImage, class TOutputImage>
const typename PreprocessingHelper<TInputImage,TOutputImage>::InputImageType *
PreprocessingHelper<TInputImage,TOutputImage>::GetQuantileThresholdedOutput() const
{
	if (!m_QuantileComputed) {
		return m_InputImage;
	}

	typename WindowFilter::Pointer f_f = WindowFilter::New();
	f_f->SetInput( m_InputImage );
	f_f->SetWindowMaximum( m_MaxThreshold );
	f_f->SetWindowMinimum( m_MinThreshold );
	f_f->SetOutputMaximum( m_Max );
	f_f->SetOutputMinimum( m_Min );
	f_f->Update();

	InputImagePointer im = f_f->GetOutput();

	return im.GetPointer();
}

template <class TInputImage, class TOutputImage>
const typename PreprocessingHelper<TInputImage,TOutputImage>::InputImageType *
PreprocessingHelper<TInputImage,TOutputImage>::GetRATSThresholdedOutput(unsigned int mode) {
	if ( m_InputBinarized.IsNotNull() ) {
		return m_InputBinarized.GetPointer();
	}

	InputPixelType insideValue = 1;

	typename BinarizeFilter::Pointer threshold = BinarizeFilter::New();
	threshold->SetInput (m_InputImage);
	threshold->SetLowerThreshold( m_RATSThreshold );
	threshold->SetInsideValue ( insideValue );
	threshold->SetOutsideValue ( 0u );
	threshold->Update();


	this->SetInputBinarized( threshold->GetOutput() );

	StrEl ball;
	typename StrEl::SizeType ballSize;
	ballSize.Fill(2.0);
	ball.SetRadius( ballSize );
	ball.CreateStructuringElement();

	typename Closing::Pointer closing = Closing::New();
	typename Opening::Pointer opening = Opening::New();
	typename Dilate::Pointer dilate = Dilate::New();

	switch( mode ) {
	case 1:
		ballSize.Fill(1.0);
		ball.SetRadius( ballSize );
		ball.CreateStructuringElement();

		opening->SetInput( threshold->GetOutput() );
		opening->SetKernel( ball );
		opening->SetForegroundValue( insideValue );
		opening->Update();

		ballSize.Fill(15.0);
		ball.SetRadius( ballSize );
		ball.CreateStructuringElement();

		closing->SetInput( opening->GetOutput() );
		closing->SetKernel( ball );
		closing->SetForegroundValue( insideValue );
		closing->Update();
		this->SetInputBinarized( closing->GetOutput() );

		break;
	case 2:
		opening->SetInput( threshold->GetOutput() );
		opening->SetKernel( ball );
		opening->SetForegroundValue( insideValue );
		opening->Update();

		closing->SetInput( opening->GetOutput() );
		closing->SetKernel( ball );
		closing->SetForegroundValue( insideValue );
		closing->Update();

		dilate->SetInput(closing->GetOutput() );
		dilate->SetForegroundValue(insideValue);
		dilate->SetKernel( ball );
		dilate->Update();

		this->SetInputBinarized( dilate->GetOutput() );

		break;

	default:

		break;
	}

	return m_InputBinarized.GetPointer();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::ApplyRATSMask(unsigned int mode) {
	if(!m_RATSComputed) {
		this->ComputeRATSThreshold(PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_POW,PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_SIGMA);
	}

	typename Masker::Pointer masker = Masker::New();
	masker->SetInput1( m_OutputImage );
	masker->SetInput2( this->GetRATSThresholdedOutput(mode) );
	masker->Update();
	m_OutputImage = masker->GetOutput();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::ApplyRATSMaskToInput(unsigned int mode) {
	if(!m_RATSComputed) {
		this->ComputeRATSThreshold(PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_POW,PreprocessingHelper<TInputImage,TOutputImage>::DEFAULT_RATS_SIGMA);
	}

	typename InputMasker::Pointer masker = InputMasker::New();
	masker->SetInput1( m_InputImage );
	masker->SetInput2( this->GetRATSThresholdedOutput(mode) );
	masker->Update();
	m_InputImage = masker->GetOutput();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SaveBinarizedOutput(const std::string file) const
{
	if (!m_QuantileComputed) {
	  return;
	}

	typename BinarizeFilter::Pointer f_f = BinarizeFilter::New();
	f_f->SetInput( m_InputImage );
	f_f->SetUpperThreshold( m_MaxThreshold );
	f_f->SetLowerThreshold( m_MinThreshold );
	f_f->SetInsideValue( 1.0 );
	f_f->SetOutsideValue( 0.0 );
	f_f->Update();

	SaveImageToFile<InputImageType> ( f_f->GetOutput(), file);
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SaveOutput(const std::string file) const
{
	SaveImageToFile<OutputImageType> ( m_OutputImage, file);
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetStripEmptyHighBins(bool collapse){
	/*
	if(!m_QuantileComputed || ?? set thresholds manually){
		// TODO
		return;
	}*/

	double binWidth = (m_Max-m_Min)/70.0;

	if ( m_MaxThreshold < (m_Max - binWidth) ) {
		InputPixelType newMax = m_MaxThreshold * 1.20;

		if (newMax > m_Max)	return;

		typename ThresholdFilter::Pointer th_filter = ThresholdFilter::New();
		th_filter->SetInput( m_OutputImage );
		th_filter->SetOutsideValue( collapse?newMax:m_Min );
		th_filter->ThresholdAbove( newMax );
		th_filter->Update();

		m_Max = newMax;
		m_OutputImage = th_filter->GetOutput();
	}
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetInputImage(const typename PreprocessingHelper::InputImageType* im){
	if (this->m_InputImage!= im) {
		this->m_InputImage= im;
		typename IOCaster::Pointer caster = IOCaster::New();
		caster->SetInput( this->m_InputImage );
		caster->Update();
		m_OutputImage = caster->GetOutput();

		typename StatsCalculator::Pointer stats = StatsCalculator::New();
		stats->SetInput( m_InputImage );
		stats->Update();

		m_Max = stats->GetMaximum();
		m_Min = stats->GetMinimum();
		m_Mean= stats->GetMean();
		m_Sigma = stats->GetSigma();

		this->Modified();
	}
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetGaussianSmoothFilter( const double variance ) {
	  typename GaussianFilter::Pointer smoother  = GaussianFilter::New();
	  smoother->SetVariance( variance );
	  smoother->SetInput( m_OutputImage );
	  smoother->Update();
	  m_OutputImage = smoother->GetOutput();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetMRIBiasFilter() {
	  typename MRIBiasFilter::Pointer bias  = MRIBiasFilter::New();
	  bias->SetInput( m_OutputImage );
	  bias->Update();
	  m_OutputImage = bias->GetOutput();
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetResample( const typename PreprocessingHelper<TInputImage,TOutputImage>::OutputImageType::SizeType size ){
	typename OutputImageType::SizeType im_size = m_OutputImage->GetLargestPossibleRegion().GetSize();
	typename OutputImageType::SpacingType im_sp = m_OutputImage->GetSpacing();
	typename OutputImageType::SpacingType new_spacing = im_sp;

	bool resizeNeeded = false;

	for( unsigned int d = 0; d< OutputImageDimension; d++){
		if ( im_size[d] > size[d] ){
			resizeNeeded = true;
			new_spacing[d] = im_sp[d]*im_size[d] / size[d];
		}
	}

	if (resizeNeeded){
		typename Resampler::Pointer r = Resampler::New();
		r->SetInput( m_OutputImage );
		r->SetOutputDirection( m_OutputImage->GetDirection() );
		r->SetOutputOrigin( m_OutputImage->GetOrigin() );
		r->SetSize( size );
		r->SetOutputSpacing( new_spacing );
		r->Update();
		m_OutputImage = r->GetOutput();
	}
}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::SetNormalization(const typename PreprocessingHelper<TInputImage,TOutputImage>::InputPixelType minValue,const typename PreprocessingHelper<TInputImage,TOutputImage>::InputPixelType maxValue){

}

template <class TInputImage, class TOutputImage>
void PreprocessingHelper<TInputImage,TOutputImage>::PrintSelf(std::ostream& os, itk::Indent indent) const {

}

#endif /* PREPROCESSINGHELPER_TXX_ */
