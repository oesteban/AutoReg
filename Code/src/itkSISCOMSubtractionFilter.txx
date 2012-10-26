/* --------------------------------------------------------------------------------------
 * File:    itkSISCOMSubtractionFilter.txx
 * Date:    02/08/2011
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2011, Oscar Esteban - oesteban@die.upm.es
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

#ifndef __itkSISCOMSubtractionFilter_txx
#define __itkSISCOMSubtractionFilter_txx

#include "itkSISCOMSubtractionFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressAccumulator.h"

#include "vnl/vnl_det.h"
#include "itkMatrix.h"


namespace itk
{

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::SISCOMSubtractionFilter(): m_UseNormalization(true)
{
	this->SetNumberOfRequiredInputs( 2 );
	//this->InPlaceOff(); // ???
	m_InputSource = this->GetInput(0);

	this->ProcessObject::SetNumberOfRequiredOutputs(1);
/*
	m_RegTransform = IdentityTransform<TInterpolatorPrecisionType, InputImageDimension>::New();
	m_InterictalTransform = IdentityTransform<TInterpolatorPrecisionType, InputImageDimension>::New();
*/
	m_IctalTransform = TransformType::New();
	m_IctalTransform->SetIdentity();

	m_Interpolator = LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>::New();
	m_DefaultPixelValue = 0;
	m_NormalizationFactor = 0.0;
	m_Resample = ResampleFilter::New();
	m_Subtract = SubtractFilter::New();
}

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::SetInterictalImage( const TInputImage * image )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<TInputImage *>( image ));
}

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::SetIctalImage( const TInputImage * image )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<TInputImage *>( image ));
}


//! Method for computing the output data
/*!
  \authors Berta Martí berta@ub.es, Oscar Esteban oesteban@die.upm.es
  \brief This method uses the Berta Martí's code in file rb9p.cpp of FocusDETLib
  \sa GetInterictalNormalizedOutput(), GetOutput()
*/
template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::GenerateData()
{
  // Resample interictal image to ictal's space
  m_Resample->SetInput( this->GetInput(0) );
  m_Resample->SetReferenceImage( this->GetInput(1) );
  m_Resample->SetUseReferenceImage(true);
  m_Resample->SetTransform( m_RegTransform );

  // These are our inputs (m_InputSource for interictal, ictalSource for ictal)
  m_InputSource = m_Resample->GetOutput();
  InputImageConstPointer ictalSource = this->GetInput(1);

  // BERTA: if I followed rb9p.cpp, here I should compute a counts factor called factorIvsII
  // but you compute it, then you apply it before registration and, you apply it inversion to the results
  // so I guess that I should not take it into accout here. If I should, here is the place.

  // Compute and apply normalization factor
  if ( m_UseNormalization ) {
	  m_Resample->Update();
	  this->ComputeNormalizationFactor();

	  typename ApplyFactorFilter::Pointer factorFilter = ApplyFactorFilter::New();
	  factorFilter->SetInput(m_InputSource);
	  factorFilter->SetConstant(m_NormalizationFactor);
	  factorFilter->Update();
	  m_InputSource = factorFilter->GetOutput();
  }

  // Get intensity maxima of the two masked images:
  InputPixelType ii_max = 0;
  InputImageConstIterator it_ii ( m_InputSource, m_InputSource->GetLargestPossibleRegion() ); it_ii.Begin();
  for ( it_ii.Begin(); !it_ii.IsAtEnd(); ++it_ii) {
	  typename InputImageType::PointType point;
	  m_InputSource->TransformIndexToPhysicalPoint( it_ii.GetIndex() , point);

  	  if ( m_IctalImageMask && !m_IctalImageMask->IsInside( point ) ) continue;

  	  if (it_ii.Get() > ii_max) ii_max=it_ii.Get();
  }

  InputPixelType i_max = 0;
  InputImageConstIterator it_i ( ictalSource, ictalSource->GetLargestPossibleRegion() );
  for ( it_i.Begin(); !it_i.IsAtEnd(); ++it_i) {
	  typename InputImageType::PointType point;
	  ictalSource->TransformIndexToPhysicalPoint( it_i.GetIndex() , point);

  	  if ( m_IctalImageMask && !m_IctalImageMask->IsInside( point ) ) continue;

  	  if (it_i.Get() > i_max) i_max=it_i.Get();
  }

  // Find uchar factor
  double factor255 = 255.0;
  factor255/= ( ii_max > i_max )?ii_max:i_max;


  std::cout <<"Norm. Factor: " << m_NormalizationFactor << ", IImax = " << ii_max << ", Imax = " << i_max << ", factor255 = " << factor255 << std::endl;


  // Apply factor to images
  // BERTA: this factor, given the fact it is the same for the 2 images, it does not affect the subtraction. IMHO, it is unnecessary
  typename ApplyFactorFilter::Pointer charFactorFilterII = ApplyFactorFilter::New();
  charFactorFilterII->SetInput( m_InputSource );
  charFactorFilterII->SetConstant( factor255 );
  charFactorFilterII->Update();

  typename ApplyFactorFilter::Pointer charFactorFilterI = ApplyFactorFilter::New();
  charFactorFilterI->SetInput( ictalSource );
  charFactorFilterI->SetConstant( factor255 );
  charFactorFilterI->Update();


  // Create and setup mean filters
  // Code copied from Berta Martí's ImaFun.cpp
  typename InputImageType::SizeType indexRadius;
  indexRadius.Fill( 1 ); // 1 pixel radius along x, y & z

  // Mean filter for interictal image
  typename MeanFilter::Pointer meanFilter1 = MeanFilter::New();
  meanFilter1->SetRadius(indexRadius);
  meanFilter1->SetInput( charFactorFilterII->GetOutput() );
  meanFilter1->Update();

  // Mean filter for ictal image
  typename MeanFilter::Pointer meanFilter2 = MeanFilter::New();
  meanFilter2->SetRadius(indexRadius);
  meanFilter2->SetInput( charFactorFilterI->GetOutput() );
  meanFilter2->Update();

  m_Subtract->SetInput1( meanFilter1->GetOutput() );
  m_Subtract->SetInput2( meanFilter2->GetOutput() );
  m_Subtract->Update();

  OutputImagePointer im_sub = m_Subtract->GetOutput();


  // Generate empty data structure for output image
  OutputImagePointer outputPtr = this->GetOutput(0);
  outputPtr->SetRegions( im_sub->GetLargestPossibleRegion() );
  outputPtr->SetSpacing( im_sub->GetSpacing() );
  outputPtr->SetDirection( im_sub->GetDirection() );
  outputPtr->SetOrigin( im_sub->GetOrigin() );
  outputPtr->Allocate();
  outputPtr->FillBuffer( 0.0 );

  // Compute II image maximum intensity
  typename CalculatorFilter::Pointer calc = CalculatorFilter::New();
  calc->SetImage( meanFilter1->GetOutput() );
  calc->Compute();
  InputPixelType max1 = calc->GetMaximum();
  InputPixelType min1 = max1 * 0.3;
  InputImagePointer im_mean1 = meanFilter1->GetOutput();


  // Iterator to copy subtraction to output
  OutputImageIterator it( outputPtr, outputPtr->GetLargestPossibleRegion());
  it.Begin();

  // Generate (copy data to output). Do not copy pixels outside masks or outside thresholds window (min1,max1)
  while ( !it.IsAtEnd() ) {
	  typename InputImageType::PointType point;
	  typename InputImageType::IndexType index = it.GetIndex();
	  ++it;
	  outputPtr->TransformIndexToPhysicalPoint( index, point);

	  if ( m_IctalImageMask && !m_IctalImageMask->IsInside( point ) ) continue;

	  typename InputImageType::PointType point_ii = m_RegTransform->TransformPoint( point );
	  if( m_InterictalImageMask && !m_InterictalImageMask->IsInside( point_ii ) ) continue;

	  typename InputImageType::IndexType idx_ii;
	  im_mean1->TransformPhysicalPointToIndex( point_ii, idx_ii );
	  InputPixelType val = im_mean1->GetPixel( idx_ii );

	  // Apply computed intensity thresholds from II
	  // BERTA: This mask seems not to affect. Also, it is not very "scientific" nor justified
	  if ( val > max1 || val < min1 ) continue;

	  it.Set( im_sub->GetPixel( it.GetIndex() ) );
  }

}

//! Method for computing the normalization factor between the two inputs
/*!
  \authors Berta Martí berta@ub.es, Oscar Esteban oesteban@die.upm.es
  \brief This method computes the normalization factor by using the Berta Martí's code
         of file fc2max.cpp
  \sa GetNormalizationFactor()
*/
template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::ComputeNormalizationFactor()
{
	if ( m_NormalizationFactor > 0.0 ) return;
	InputImageConstPointer ictalImage = this->GetInput(1);

	//Inital variables
	int num;
	double nint = 100.0;

	//Histogram length
	int NINT=100;
	int P=10*NINT-1;

	typedef std::vector <double> VectorDouble;
	//Initialize the vectors of histogram (l and h)
	VectorDouble l;
	l.resize(P);
	VectorDouble h;
	h.resize(P);

	//Profile ratios to get the inital factor (I/II)
	typedef itk::ImageRegionConstIterator<InputImageType> Iterator;
	Iterator it(ictalImage, ictalImage->GetLargestPossibleRegion());

	for (it = it.Begin(); !it.IsAtEnd(); ++it) {
		typename InputImageType::IndexType idx = it.GetIndex();
		typename InputImageType::PointType point;
		ictalImage->TransformIndexToPhysicalPoint( idx, point);

		if ( m_IctalImageMask.IsNotNull() && !m_IctalImageMask->IsInside( point ) ) continue;

		InputPixelType pixelII = m_InputSource->GetPixel(idx);
		if ( pixelII <= 0 ) continue;
		double ratio = ( ictalImage->GetPixel(idx) / pixelII );

		if (ratio>0.) {
			num = (int) floor( ratio*nint+.5);
			if (num<P) h[num]+=1;
		}
	}

	//Smooth histogram
	for (int i=2;i<P-2;i++) {
		l[i]=(h[i-2]+h[i-1]+h[i]+h[i+1]+h[i+2])/5.;
	}

	// Get the maximum of l vector and its position
	double maxValue = *std::max_element(l.begin(), l.end());
	VectorDouble::iterator iterator = max_element(l.begin(), l.end());
	int positionMax = std::distance(l.begin(), iterator);

	//Average
	double N=0.0;
	double sx=0.0;
	for (int i=positionMax-NINT/20; i<=positionMax+NINT/20; i++) {
		N += h[i];
		sx += i*h[i]/nint;
	}

	double average = sx/N;


	//Fitting a parabola  y=ax2+bx+c
	//Init variables
	double sx1=0.0;
	double sx2=0.0;
	double sx3=0.0;
	double sx4=0.0;
	double sy=0.0;
	double sxy=0.0;
	double sx2y=0.0;

	double ix;
	for(int i=positionMax-NINT/10;i<positionMax+NINT/10;i++){
		ix=(float)(i);
		sx1+=ix;
		sx2+=ix*ix;
		sx3+=ix*ix*ix;
		sx4+=ix*ix*ix*ix;
		sy+=h[i];
		sxy+=ix*h[i];
		sx2y+=ix*ix*h[i];
	}

	double n_p=nint/5.0;

	typedef itk::Matrix< double, 3, 3 > Matrix;
	Matrix S; S.Fill( 0.0 );
	S(0,0) = sx4; S(0,1) = sx3; S(0,2) = sx2;
	S(1,0) = sx3; S(1,1) = sx2; S(1,2) = sx1;
	S(2,0) = sx2; S(2,1) = sx1; S(2,2) = n_p;
	double denominador = vnl_det( S.GetVnlMatrix());

	Matrix A(S); A(0,0) = sx2y; A(1,0) = sxy; A(2,0) = sy;
	double a = vnl_det( A.GetVnlMatrix() );

	Matrix B(S); B(0,1) = sx2y; B(1,1) = sxy; B(2,1) = sy;
	double b = vnl_det( B.GetVnlMatrix() );

	Matrix C(S); C(0,2) = sx2y; C(1,2) = sxy; C(2,2) = sy;
	double c = vnl_det( C.GetVnlMatrix() );

	a/=denominador;
	b/=denominador;
	c/=denominador;

	m_NormalizationFactor =-(b)/(2.*a)/nint;

	double new_h=-b*b/(4.*a)+c;
}

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::ComposeTransforms()
{
	typename InputImageType::PointType o_bas = m_InputSource->GetOrigin();
	typename InputImageType::PointType o_act = this->GetInput(1)->GetOrigin();

	typename InputImageType::PointType o_zero; o_zero.Fill(0.0);

	typename TransformType::MatrixType M_reg_inv = m_RegTransform->GetInverseMatrix();
	typename TransformType::OffsetType O_reg = m_RegTransform->GetOffset();

	typename TransformType::MatrixType M_coreg = m_InterictalTransform->GetMatrix();
	typename TransformType::OffsetType O_coreg = m_InterictalTransform->GetOffset();

	typename TransformType::MatrixType M_final = M_reg_inv * M_coreg;
	typename TransformType::OffsetType O_final = M_reg_inv * (O_coreg - O_reg - (o_bas - o_zero)) + (o_act - o_zero);

	m_IctalTransform = TransformType::New();
	m_IctalTransform->SetMatrix(M_final);
	m_IctalTransform->SetOffset(O_final);
}

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::ResampleToReferenceSpace( const typename SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>::InputImageType* i_ref )
 {
	this->ComposeTransforms();

	typename ResampleFilter::Pointer res_bas = ResampleFilter::New();
	res_bas->SetUseReferenceImage(true);
	res_bas->SetReferenceImage(i_ref);
	res_bas->SetInput( this->GetInterictalNormalizedOutput() );
	res_bas->SetTransform(m_IctalTransform);
	res_bas->Update();
	m_InterictalRefSpace = res_bas->GetOutput();

	typename ResampleFilter::Pointer res_act = ResampleFilter::New();
	res_act->SetUseReferenceImage(true);
	res_act->SetReferenceImage(i_ref);
	res_act->SetInput( this->GetInput(1) );
	res_act->SetTransform(m_IctalTransform);
	res_act->Update();
	m_IctalRefSpace = res_act->GetOutput();

	typedef itk::ResampleImageFilter< OutputImageType, OutputImageType > DiffResampler;
	typename DiffResampler::Pointer res_dif = DiffResampler::New();
	res_dif->SetOutputParametersFromImage( i_ref );
	res_dif->SetInput( this->GetOutput() );
	res_dif->SetTransform( m_IctalTransform );
	res_dif->Update();
	m_SubtractionRefSpace = res_dif->GetOutput();
 }

template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
typename SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>::TransformPointerType
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::GetIctalTransform() {
	this->ComposeTransforms();
	return m_IctalTransform;
}


template< class TInputImage, class TOutputImage, class TInterpolatorPrecisionType >
void
SISCOMSubtractionFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "UseNormalization?: " << m_UseNormalization << std::endl;
  os << indent << "Registration Transform: " << m_RegTransform.GetPointer() << std::endl;
  os << indent << "Co-registration Transform: " << m_InterictalTransform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "DefaultPixelValue: "
     << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_DefaultPixelValue)
     << std::endl;
}

} // end namespace itk


#endif
