/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToAnisotropicHistogramGenerator.txx,v $
  Language:  C++
  Date:      $Date: 2009-08-17 18:29:01 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarImageToAnisotropicHistogramGenerator_txx
#define __itkScalarImageToAnisotropicHistogramGenerator_txx

#include "itkScalarImageToAnisotropicHistogramGenerator.h"


namespace itk { 
namespace Statistics {


template < class TImage >
ScalarImageToAnisotropicHistogramGenerator< TImage >
::ScalarImageToAnisotropicHistogramGenerator() 
{
  m_ImageToListAdaptor = AdaptorType::New();
  m_HistogramGenerator = GeneratorType::New();
  m_HistogramGenerator->SetInput( m_ImageToListAdaptor );
}

template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetInput( const ImageType * image ) 
{
  m_ImageToListAdaptor->SetImage( image );
}


template < class TImage >
const typename ScalarImageToAnisotropicHistogramGenerator< TImage >::HistogramType *
ScalarImageToAnisotropicHistogramGenerator< TImage >
::GetOutput() const
{
  return m_HistogramGenerator->GetOutput();
}

template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::Compute() 
{
  m_HistogramGenerator->Update();
}

template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetNumberOfBins( unsigned int numberOfBins ) 
{
  typename HistogramType::SizeType size;
  size.SetSize(1);
  size.Fill( numberOfBins );
  m_HistogramGenerator->SetHistogramSize( size );
}


template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetHistogramMin( RealPixelType minimumValue ) 
{
  typedef typename GeneratorType::HistogramMeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType minVector;
  minVector[0] = minimumValue;
  m_HistogramGenerator->SetHistogramBinMinimum( minVector );
}


template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetMaximumNumberOfIterations( unsigned int n )
{
  m_HistogramGenerator->SetMaxEstimatorIterations( n );
}


template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetHistogramMax( RealPixelType maximumValue ) 
{
  typedef typename GeneratorType::HistogramMeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType maxVector;
  maxVector[0] = maximumValue;
  m_HistogramGenerator->SetHistogramBinMaximum( maxVector );
}

template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::SetMarginalScale( double marginalScale )
{
  m_HistogramGenerator->SetMarginalScale( marginalScale );
}

template < class TImage >
void
ScalarImageToAnisotropicHistogramGenerator< TImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << "ImageToListSample adaptor = " << m_ImageToListAdaptor << std::endl;
  os << "HistogramGenerator = " << m_HistogramGenerator << std::endl;
}

} // end of namespace Statistics 
} // end of namespace itk

#endif
