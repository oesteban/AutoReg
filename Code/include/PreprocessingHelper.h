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

#ifndef PREPROCESSINGHELPER_H_
#define PREPROCESSINGHELPER_H_

#include <itkObject.h>
#include <itkSmartPointer.h>

#include <itkNumericTraits.h>

#include <itkCastImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkOtsuThresholdImageCalculator.h>
#include <itkAdaptativeBiningThresholdImageCalculator.h>
#include <itkQuantileThresholdImageCalculator.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkRobustAutomaticThresholdCalculator.h>

#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkMRIBiasFieldCorrectionFilter.h>
#include <itkHistogramMatchingImageFilter.h>

template <class TInputImage, class TOutputImage>
class PreprocessingHelper : public itk::Object
{

  public:
    typedef PreprocessingHelper                         Self;
    typedef itk::Object                                 Superclass;
    typedef itk::SmartPointer< Self >                   Pointer;
    typedef itk::SmartPointer< const Self >             ConstPointer;

    itkTypeMacro( PreprocessingHelper, Object );
    itkNewMacro( Self );

    typedef TInputImage                                 InputImageType;
    typedef typename InputImageType::Pointer            InputImagePointer;
    typedef typename InputImageType::ConstPointer       InputImageConstPointer;
    typedef typename InputImageType::PixelType          InputPixelType;

    typedef TOutputImage                                OutputImageType;
    typedef typename OutputImageType::Pointer           OutputImagePointer;
    typedef typename OutputImageType::ConstPointer      OutputImageConstPointer;
    typedef typename OutputImageType::PixelType         OutputPixelType;

    itkStaticConstMacro(InputImageDimension, unsigned int, InputImageType::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int, OutputImageType::ImageDimension);

    const static float DEFAULT_RATS_POW;
    const static float DEFAULT_RATS_SIGMA;

    typedef typename itk::StatisticsImageFilter<InputImageType>   StatsCalculator;
    typedef typename itk::StatisticsImageFilter<OutputImageType>  OutCalculator;

    typedef typename itk::BinaryThresholdImageFilter
                            <InputImageType, InputImageType>      BinarizeFilter;

	typedef typename itk::BinaryBallStructuringElement< InputPixelType, OutputImageDimension>   StrEl;
	typedef typename itk::BinaryMorphologicalOpeningImageFilter< InputImageType, InputImageType, StrEl >  Opening;
	typedef typename itk::BinaryMorphologicalClosingImageFilter< InputImageType, InputImageType, StrEl >  Closing;
	typedef typename itk::BinaryDilateImageFilter< InputImageType, InputImageType, StrEl > Dilate;

	typedef typename itk::IntensityWindowingImageFilter
	                    <InputImageType,InputImageType>           WindowFilter;
	typedef typename itk::QuantileThresholdImageCalculator
	                                           <InputImageType>   QuantileCalculator;

	typedef typename itk::GradientMagnitudeRecursiveGaussianImageFilter
	                         < InputImageType, InputImageType >   GradientType;
	typedef typename GradientType::OutputImageType                GradientOutputType;

	typedef typename itk::RobustAutomaticThresholdCalculator
	                      < InputImageType, GradientOutputType>   CalculatorType;

    typedef typename itk::CastImageFilter
                        < InputImageType, OutputImageType>        IOCaster;

    typedef typename itk::DiscreteGaussianImageFilter
                        < OutputImageType, OutputImageType >      GaussianFilter;

    typedef typename itk::MRIBiasFieldCorrectionFilter
             < OutputImageType, OutputImageType, InputImageType > MRIBiasFilter;

    typedef typename itk::ResampleImageFilter
                        < OutputImageType, OutputImageType >      Resampler;

	typedef typename itk::ThresholdImageFilter<OutputImageType>   ThresholdFilter;


	typedef typename itk::MaskImageFilter< OutputImageType ,
                               InputImageType, OutputImageType >   Masker;
	typedef typename itk::MaskImageFilter< InputImageType ,
                               InputImageType, InputImageType >    InputMasker;

	typedef typename itk::HistogramMatchingImageFilter
	           < OutputImageType,OutputImageType >				   HistogramFilter;


    /** Get the output image */
    itkGetConstObjectMacro(OutputImage,OutputImageType );

    itkGetConstObjectMacro(InputBinarized,InputImageType );
    itkSetConstObjectMacro(InputBinarized,InputImageType );

    itkSetClampMacro( MinThreshold, InputPixelType,
    		itk::NumericTraits<InputPixelType>::min(),
    		itk::NumericTraits<InputPixelType>::max() );

	itkGetConstMacro( MinThreshold, InputPixelType );

    itkSetClampMacro( MaxThreshold, InputPixelType,
    		itk::NumericTraits<InputPixelType>::min(),
    		itk::NumericTraits<InputPixelType>::max() );

	itkGetConstMacro( MaxThreshold, InputPixelType );

    itkSetClampMacro( RATSThreshold, InputPixelType,
    		itk::NumericTraits<InputPixelType>::min(),
    		itk::NumericTraits<InputPixelType>::max() );

	itkGetConstMacro( RATSThreshold, InputPixelType );

	itkGetConstMacro( Min, InputPixelType);
	itkGetConstMacro( Max, InputPixelType);
	itkGetConstMacro( Mean, InputPixelType);
	itkGetConstMacro( Sigma, InputPixelType);
	itkGetConstMacro( Offset, InputPixelType);

	void ComputeQuantileThresholds(float minPercentage, float maxPercentage );
	void ComputeRATSThreshold(float pow, float sigma );
	void ApplyRATSMask(unsigned int mode);
	void ApplyRATSMaskToInput(unsigned int mode);
	void SetHistogramMatchingImage(const OutputImageType * refImage);

	virtual const InputImageType * GetQuantileThresholdedOutput() const;
	virtual const InputImageType * GetRATSThresholdedOutput(unsigned int mode);

	void SaveBinarizedOutput(const std::string file) const;
	void SaveOutput(const std::string file) const;

	virtual void SetInputImage(const InputImageType* im);

	void SetGaussianSmoothFilter( const double variance );
	void SetMRIBiasFilter( );

	void SetResample( const typename OutputImageType::SizeType size );

	void SetStripEmptyHighBins( bool collapse=true );

	void SetNormalization(const InputPixelType minValue,const InputPixelType maxValue);

	// void SetIntensityNormalization( OutputImageType::PixelType min, OutputImageType::PixelType max );

  protected:
    PreprocessingHelper();
    ~PreprocessingHelper() {};
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    PreprocessingHelper(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputImageConstPointer    m_InputImage;
    OutputImageConstPointer   m_OutputImage;
    InputImageConstPointer    m_InputBinarized;

    InputPixelType m_MaxThreshold;
    InputPixelType m_MinThreshold;
    InputPixelType m_Max;
    InputPixelType m_Min;
    InputPixelType m_Mean;
    InputPixelType m_Sigma;
    InputPixelType m_Offset;
    InputPixelType m_RATSThreshold;

    bool                      m_QuantileComputed;
    bool					  m_RATSComputed;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "PreprocessingHelper.txx"
#endif

#endif /* PREPROCESSINGHELPER_H_ */
