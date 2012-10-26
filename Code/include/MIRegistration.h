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

#ifndef MIREGISTRATION_H
#define MIREGISTRATION_H

#ifdef USE_ADAPTATIVE_BINING
  #include "../Review/itkHistogramImageToImageMetric.h"
#endif

#include <itkObject.h>
#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkImageRegistrationMethod.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkCastImageFilter.h>

// PrepareRegImage includes start
#include <itkMinimumMaximumImageCalculator.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkScalarImageToHistogramGenerator.h>
// PrepareRegImage includes end


template< class TImage>
static void SaveImageToFile( typename TImage::ConstPointer im, std::string name )
{
  typedef itk::ImageFileWriter< TImage >  Writer;
  typename Writer::Pointer w = Writer::New();
  w->SetInput( im );
  w->SetFileName( name.c_str() );
  w->Update();
}

#include "PreprocessingHelper.h"
#include "OptimizerHelper.h"
#include "TransformHelper.h"
#include "MetricHelper.h"



template <class TFixedImage, class TMovingImage, class TScalarType=double >
class MIRegistration: public itk::Object
{
  
  public:
    
    typedef MIRegistration                         Self;
    typedef Object                                 Superclass;
    typedef itk::SmartPointer< Self >              Pointer;
    typedef itk::SmartPointer< const Self >        ConstPointer;

    itkTypeMacro( MIRegistration, Object );

    itkNewMacro( Self );
    
    /** Input Images definitions */
    typedef TFixedImage                                 InputFixedImageType;
    typedef typename InputFixedImageType::Pointer       InputFixedImagePointer;
    typedef typename InputFixedImageType::ConstPointer  InputFixedImageConstPointer;
    typedef typename InputFixedImageType::RegionType    InputFixedImageRegionType;
    typedef typename InputFixedImageType::PixelType     InputFixedImagePixelType;
    typedef TMovingImage                                InputMovingImageType;
    typedef typename InputMovingImageType::Pointer      InputMovingImagePointer;
    typedef typename InputMovingImageType::ConstPointer InputMovingImageConstPointer;
    typedef typename InputMovingImageType::RegionType   InputMovingImageRegionType;
    typedef typename InputMovingImageType::PixelType    InputMovingImagePixelType;

    /** Internal images definitions */
    itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);
    itkStaticConstMacro(MovingImageDimension, unsigned int, TMovingImage::ImageDimension);
    typedef TScalarType                                                           InternalPixelType;
    typedef typename itk::Image
            <InternalPixelType, itkGetStaticConstMacro(FixedImageDimension)>      InternalFixedImageType;
    typedef typename itk::Image
            <InternalPixelType, itkGetStaticConstMacro(MovingImageDimension)>     InternalMovingImageType;
    typedef typename InternalFixedImageType::Pointer                              InternalFixedPointer;
    typedef typename InternalFixedImageType::ConstPointer                         InternalFixedConstPointer;
    typedef std::vector< InternalFixedConstPointer >                              InternalFixedImageList;
    typedef typename InternalMovingImageType::Pointer                             InternalMovingPointer;
    typedef typename InternalMovingImageType::ConstPointer                        InternalMovingConstPointer;
    typedef std::vector< InternalMovingConstPointer >                             InternalMovingImageList;

    typedef typename itk::CastImageFilter<InputFixedImageType,InternalFixedImageType> InternalFixedCaster;
    typedef typename itk::CastImageFilter<InputMovingImageType,InternalMovingImageType> InternalMovingCaster;
    
    typedef typename itk::MaskImageFilter< InternalFixedImageType, InputFixedImageType, InternalFixedImageType >        InternalFixedImageMasker;
    typedef typename itk::MaskImageFilter< InternalMovingImageType, InputMovingImageType, InternalMovingImageType >     InternalMovingImageMasker;

    typedef itk::LinearInterpolateImageFunction<InternalFixedImageType, double>   Interpolator;
    typedef itk::ResampleImageFilter
                           < InputMovingImageType, InputFixedImageType >          ResampleFilter;

    typedef itk::ImageRegistrationMethod
                        < InternalFixedImageType, InternalMovingImageType >       RegistrationType;
    typedef typename RegistrationType::Pointer                                    RegistrationPointer;
    typedef std::vector< RegistrationPointer >                                    RegistrationPointerList;
				  
    typedef itk::MultiResolutionPyramidImageFilter
                        < InternalFixedImageType, InternalMovingImageType >       PyramidFilter;

    typedef PreprocessingHelper< InputFixedImageType, InternalFixedImageType>     FixedPHelper;
    typedef typename FixedPHelper::Pointer                                        FixedPHelperPointer;
    typedef std::vector< FixedPHelperPointer >                                    FixedPHelperList;

    typedef PreprocessingHelper< InputMovingImageType, InternalMovingImageType>   MovingPHelper;
    typedef typename MovingPHelper::Pointer                                       MovingPHelperPointer;
    typedef std::vector< MovingPHelperPointer >                                   MovingPHelperList;

    typedef OptimizerHelper                                                       OHelper;
    typedef OptimizerHelper::Pointer                                              OptimizerHelperPointer;
    typedef std::vector< OptimizerHelperPointer >                                 OptimizerHelperList;
    typedef OptimizerHelper::OPTIMIZER_TYPE                                       OptimizerTypeEnum;
    typedef typename OptimizerHelper::GenericType::Pointer                        OptimizerPointer;
    typedef OHelper::ScalesType                                                   ScalesType;

    typedef TransformHelper
          < InternalFixedImageType, InternalMovingImageType, double >             THelper;
    typedef typename THelper::Pointer                                             TransformHelperPointer;
    typedef std::vector< TransformHelperPointer >                                 TransformHelperList;
    typedef typename THelper::TRANSFORM_TYPE                                      TransformTypeEnum;
    typedef typename THelper::INITIALIZATION_TYPE                                 InitializationTypeEnum;
//     typedef typename THelper::GenericType::Pointer                                TransformPointer;
    typedef typename THelper::GenericType*                                        TransformPointer;

    typedef MetricHelper< InternalFixedImageType, InternalMovingImageType>        MHelper;
    typedef typename MHelper::Pointer                                             MetricHelperPointer;
    typedef std::vector< MetricHelperPointer >                                    MetricHelperList;
    typedef typename MHelper::METRIC_TYPE                                         MetricTypeEnum;
    typedef typename MHelper::GenericType::Pointer                                MetricPointer;

    itkSetMacro( OutputDirectory, std::string );
    itkGetConstReferenceMacro( OutputDirectory, std::string);
    
    itkSetMacro( SaveRegImages, bool );
    itkGetConstReferenceMacro( SaveRegImages, bool );
    
    itkSetMacro(InitFixedLowerBound, typename InputFixedImageType::PixelType );
    itkGetConstReferenceMacro( InitFixedLowerBound, typename InputFixedImageType::PixelType );

    itkSetMacro(InitMovingLowerBound, typename InputMovingImageType::PixelType );
    itkGetConstReferenceMacro( InitMovingLowerBound, typename  InputMovingImageType::PixelType );
    
    itkSetMacro( UseBinarizedInit, bool );
    itkGetConstReferenceMacro( UseBinarizedInit, bool );

    itkSetMacro( UseFixedOriginalInit, bool );
    itkGetConstReferenceMacro( UseFixedOriginalInit, bool );

    itkSetMacro( UseMovingOriginalInit, bool );
    itkGetConstReferenceMacro( UseMovingOriginalInit, bool );

    itkSetMacro( InitSize, unsigned int );
    itkGetConstMacro( InitSize, unsigned int);

    FixedPHelperPointer GetFixedPreprocessHelper( const unsigned int level) const
    { return m_FixedPreprocessHelpers[level];   }
    MovingPHelperPointer GetMovingPreprocessHelper( const unsigned int level) const
    { return m_MovingPreprocessHelpers[level];   }
    
    void SetOutputStream( std::ostream* os );
    
    void SetFixedImage( InputFixedImageConstPointer fixedImage );

    void SetMovingImage( InputMovingImageConstPointer movingImage );
		       
    void NormalizeImages( unsigned int level, double minValue, double maxValue );
    void SmoothImages( unsigned int level, double smoothFixed, double smoothMoving );
    void ResampleImages( unsigned int level, unsigned int fixedNOPD, unsigned int movingNOPD );
    
    void SetFixedMovingSwapForLevel( unsigned int level, bool value= true);

    void SetMultiResolutionLevels( unsigned int levels );

    void SetFixedImageMask( InputFixedImageConstPointer fixedMask )
    { m_FixedMask = fixedMask; }
    
    void SetMovingImageMask( InputMovingImageConstPointer movingMask )
    { m_MovingMask = movingMask; }

    void SetOptimizerType( OptimizerTypeEnum type, unsigned int level );
    void SetOptimizerType( OptimizerTypeEnum type );
    OptimizerHelperPointer GetOptimizerHelper( unsigned int level = 0 )
    { return m_OptimizerHelpers[level]; }
    
    void SetOptimizerParameters( unsigned int level, OptimizerHelper::ParametersType param );
    void SetOptimizerHelperParametersUpdate( unsigned int level, OptimizerHelper::ParametersType param );
    
    MetricHelperPointer GetMetricHelper( unsigned int level )
    { return m_MetricHelpers[level]; }
    
    
    void SetMetricParameters( unsigned int level, typename MHelper::ParametersType param );
    void SetMetricParametersUpdate( unsigned int level, typename MHelper::ParametersType param );
    void SetMetricUseMasks( unsigned int level, bool useMasks = true)
    { m_UseMasksRegistrationLevel[level] = useMasks; }
    
    void SetInitializationType( InitializationTypeEnum type );
    void SetInitialTransform( TransformPointer transform );
    void SetInitializationUseMasks( bool useMasks = true)
    { m_UseMasksInitialization = useMasks; }
    
    /*void SetInitializationUseBinarizedImages( bool useBinarized = true );*/
    
    void SetTransformType( TransformTypeEnum type, unsigned int level );
    void SetTransformType( TransformTypeEnum type );
    TransformHelperPointer GetTransformHelper( unsigned int level = 0 )
    { return m_TransformHelpers[level]; }
    
    int StartRegistration();
    void SaveResampledImage();
    
    template <typename TImage>
    ScalesType GetAutoScales(const TImage* image, bool normalize = true) const;
    
  protected:
    MIRegistration();
    ~MIRegistration() {};
    
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
    
  private:
    unsigned int                    m_Levels;
    std::vector <bool>              m_SwapLevels;
    bool                            m_UseMasksInitialization;
    std::vector <bool>              m_UseMasksRegistrationLevel;
    bool                            m_UseBinarizedInit;
    bool                            m_UseFixedOriginalInit;
    bool                            m_UseMovingOriginalInit;
    unsigned int					m_InitSize;
    
    InternalFixedImageList          m_FixedImage;
    InternalMovingImageList         m_MovingImage;
    
    InputMovingImageConstPointer    m_InputMovingImage;
    InputFixedImageConstPointer     m_InputFixedImage;
    
    InputMovingImageConstPointer    m_MovingMask;
    InputFixedImageConstPointer     m_FixedMask;
    
    FixedPHelperList                m_FixedPreprocessHelpers;
    MovingPHelperList               m_MovingPreprocessHelpers;
    OptimizerHelperList             m_OptimizerHelpers;
    MetricHelperList                m_MetricHelpers;
    TransformHelperList             m_TransformHelpers;
    RegistrationPointerList         m_RegistrationPointers;
    
    std::string                     m_OutputDirectory;
    bool                            m_SaveRegImages;
    std::ostream*                   m_OutputStream;
    
    typename InputFixedImageType::PixelType   m_InitFixedLowerBound;
    typename InputMovingImageType::PixelType  m_InitMovingLowerBound;
    
    void Connect(unsigned int level);
    
    template <class TInputImage, class TOutputImage >
    typename TOutputImage::Pointer PrepareRegImage (typename TInputImage::ConstPointer image,
						    float thPercentage,
						    double rescaleMax );
						    
    template <class TImage>
    typename TImage::Pointer NormalizeImage ( typename TImage::ConstPointer image, typename TImage::PixelType minValue, typename TImage::PixelType maxValue );
    
    template <class TImage>
    typename TImage::Pointer SmoothImage ( typename TImage::ConstPointer image, double smooth );
    
    template <class TImage>
    typename TImage::Pointer ResampleImage ( typename TImage::ConstPointer image, unsigned int numberOfPixelsOfDimension );
    
    void InitImagesProcessing(bool binarize=true);
};


#include "MIRegistration.txx"


#endif // MIREGISTRATION_H
