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

#ifndef METRICHELPER_H
#define METRICHELPER_H

#include <boost/any.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <stdio.h>
#include <iostream>


#include <itkObject.h>
#include <itkImageMaskSpatialObject.h>

#ifdef USE_ADAPTATIVE_BINING
  #include "../Review/itkHistogramImageToImageMetric.h"
#endif

#include <itkArray.h>
#include <itkImageToImageMetric.h>
#include <itkMattesMutualInformationImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkMutualInformationHistogramImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>

template <class TFixedImage, class TMovingImage>
class MetricHelper : public itk::Object
{
  public:
    typedef MetricHelper                                Self;
    typedef itk::Object                                 Superclass;
    typedef itk::SmartPointer< Self >                   Pointer;
    typedef itk::SmartPointer< const Self >             ConstPointer;

    itkTypeMacro( MetricHelper, Object );

    itkNewMacro( Self );
    
    /** Some convenient typedefs. */
    typedef TFixedImage                            FixedImageType;
    typedef typename FixedImageType::Pointer       FixedImagePointer;
    typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;
    typedef typename FixedImageType::PixelType	   FixedPixelType;
    typedef TMovingImage                           MovingImageType;
    typedef typename MovingImageType::Pointer      MovingImagePointer;
    typedef typename MovingImageType::ConstPointer MovingImageConstPointer;
    typedef typename MovingImageType::PixelType	   MovingPixelType;

    itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);
    itkStaticConstMacro(MovingImageDimension, unsigned int, TMovingImage::ImageDimension);
    

    typedef typename itk::ImageMaskSpatialObject< FixedImageDimension >			  FixedMaskSpatialObjectType;
    typedef typename FixedMaskSpatialObjectType::Pointer					      FixedMaskPointer;
    typedef typename FixedMaskSpatialObjectType::ConstPointer					  FixedMaskConstPointer;
    typedef typename FixedMaskSpatialObjectType::ImageType						  FixedMaskImageType;
    typedef const FixedMaskImageType*											  FixedMaskImagePointer;

    typedef typename itk::ImageMaskSpatialObject< MovingImageDimension >		  MovingMaskSpatialObjectType;
    typedef typename MovingMaskSpatialObjectType::Pointer						  MovingMaskPointer;
    typedef typename MovingMaskSpatialObjectType::ConstPointer					  MovingMaskConstPointer;
    typedef typename MovingMaskSpatialObjectType::ImageType					  	  MovingMaskImageType;
    typedef const MovingMaskImageType*												  MovingMaskImagePointer;

    typedef typename itk::CastImageFilter<FixedImageType, FixedMaskImageType >    FixedMaskCasterType;
    typedef typename itk::CastImageFilter<MovingImageType, MovingMaskImageType >  MovingMaskCasterType;

    /** Get the size of the input space */
    unsigned int GetFixedImageDimension(void) const {return FixedImageDimension;}

    /** Get the size of the output space */
    unsigned int GetMovingImageDimension(void) const {return MovingImageDimension;}

    enum METRIC_TYPE {
      MI,
      NMI,
      MATTES,
      VIOLA_WELLS
    };
    typedef itk::ImageToImageMetric
                     < FixedImageType, MovingImageType>                    GenericType;
    typedef itk::MattesMutualInformationImageToImageMetric
                     < FixedImageType, MovingImageType>                    MattesType;
    typedef itk::MutualInformationImageToImageMetric
                     < FixedImageType, MovingImageType>                    ViolaWellsType;
    typedef itk::MutualInformationHistogramImageToImageMetric
                     < FixedImageType, MovingImageType>                    MIType;
    typedef itk::NormalizedMutualInformationHistogramImageToImageMetric
                     < FixedImageType, MovingImageType>                    NMIType;

    typedef std::vector< boost::any >                                      ParametersType;
    typedef itk::Array< double >                                           ScalesType;
    
    
    template <class TRegistration>
    void Connect( TRegistration* r);

    void SetMetricType( METRIC_TYPE type );
    
    void SetDerivativeStepLengthScales( ScalesType scales )
    { m_Scales = scales; }
    
    void SetDerivativeStepLength( double arg )
    { m_StepLength = arg; }
    
    template < class TMetric >
    TMetric* GetMetric();
    
    void SetParameters( ParametersType param )
    { m_Parameters = param; }
    
    unsigned int GetNumberOfParameters();
    
    ParametersType GetParameters()
    { return m_Parameters; }
    
    void SetFixedImageMask( const FixedImageType* fixedMask );
    void SetMovingImageMask( const MovingImageType* movingMask );

  protected:
    MetricHelper();
    ~MetricHelper() {};
    
   virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
    
  private:
    METRIC_TYPE                          m_Type;
    ParametersType                       m_Parameters;
    ScalesType                           m_Scales;
    double                               m_StepLength;
/*
    FixedImagePointer                    m_FixedImage;
    MovingImagePointer                   m_MovingImage;
*/

    FixedMaskPointer                m_FixedMask;
    MovingMaskPointer               m_MovingMask;

    typename GenericType::Pointer        m_Metric;
    typename MattesType::Pointer         m_MattesMetric;
    typename ViolaWellsType::Pointer     m_ViolaWellsMetric;
    typename MIType::Pointer             m_MIMetric;
    typename NMIType::Pointer            m_NMIMetric;
};

#include "MetricHelper.txx"

#endif // METRICHELPER_H
