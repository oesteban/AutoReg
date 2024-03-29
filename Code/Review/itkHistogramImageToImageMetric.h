/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHistogramImageToImageMetric.h,v $
  Language:  C++
  Date:      $Date: 2009-05-02 05:43:54 $
  Version:   $Revision: 1.17 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkHistogramImageToImageMetric_h

#ifndef _WIN32
	#warning "Using ReviewHistogramImageToImageMetric"
#endif

#define __itkHistogramImageToImageMetric_h

#include "itkHistogram.h"
#include "itkImageToImageMetric.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkSampleToAnisotropicHistogramFilter.h"

namespace itk
{
/** \class HistogramImageToImageMetric 
    \brief Computes similarity between two objects to be registered
 
  This class is templated over the type of the fixed and moving
  images to be compared.
 
  The metric computes the similarity measure between pixels in the
  moving image and pixels in the fixed image using a histogram.
 
  \ingroup RegistrationMetrics */
template <class TFixedImage, class TMovingImage>
class ITK_EXPORT HistogramImageToImageMetric : 
public ImageToImageMetric<TFixedImage, TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef HistogramImageToImageMetric                   Self;
  typedef ImageToImageMetric<TFixedImage, TMovingImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(HistogramImageToImageMetric, ImageToImageMetric);
 
  /** Types transferred from the base class */
  typedef typename Superclass::RealType                   RealType;
  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::TransformParametersType
                                                          TransformParametersType;
  typedef typename Superclass::TransformJacobianType
                                                          TransformJacobianType;
  typedef typename Superclass::GradientPixelType          GradientPixelType;
  typedef typename Superclass::InputPointType             InputPointType;
  typedef typename Superclass::OutputPointType            OutputPointType;
  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::FixedImageType             FixedImageType;
  typedef typename Superclass::FixedImageType::PixelType  FixedImagePixelType;
  typedef typename Superclass::MovingImageType            MovingImageType;
  typedef typename Superclass::MovingImageType::PixelType MovingImagePixelType;
  typedef typename Superclass::FixedImageConstPointer
                                                          FixedImageConstPointerType;
  typedef typename Superclass::MovingImageConstPointer
                                                          MovingImageConstPointerType;

  /** Typedefs for histogram. This should have been defined as
      Histogram<RealType,2> but a bug in VC++7 produced an internal compiler
      error with such declaration. */
#ifdef ITK_USE_REVIEW_STATISTICS
  typedef Statistics::Histogram<double>                  HistogramType;
#else
  typedef Statistics::Histogram<double, 2>               HistogramType;
#endif

  typedef typename HistogramType::MeasurementVectorType  MeasurementVectorType;
  typedef typename HistogramType::SizeType               HistogramSizeType;
  typedef typename HistogramType::Pointer                HistogramPointer;
  typedef typename HistogramType::BinMinContainerType    BinMinContainerType;
  typedef typename HistogramType::BinMaxContainerType    BinMaxContainerType;
  

  
  typedef itk::Statistics::ImageToListSampleAdaptor<TFixedImage>     FixedAdaptorType;
  typedef itk::Statistics::ImageToListSampleAdaptor<TMovingImage>    MovingAdaptorType;
  typedef itk::Statistics::SampleToAnisotropicHistogramFilter< FixedAdaptorType, HistogramType >    FixedHistogramFilter;
  typedef itk::Statistics::SampleToAnisotropicHistogramFilter< MovingAdaptorType, HistogramType >    MovingHistogramFilter;


  /** Initializes the metric. */
  void Initialize() throw (ExceptionObject);

  /** Define the transform and thereby the parameter space of the metric
   *   and the space of its derivatives */
  void SetTransform( TransformType * transform );

  /** Sets the histogram size. Note this function must be called before
      \c Initialize(). */
  itkSetMacro( HistogramSize, HistogramSizeType );

  /** Gets the histogram size. */
  itkGetConstReferenceMacro( HistogramSize, HistogramSizeType );

  /** Factor to increase the upper bound for the samples in the histogram.
      Default value is 0.001 */
  itkSetMacro( UpperBoundIncreaseFactor, double );
  itkGetConstMacro( UpperBoundIncreaseFactor, double );

  /** The padding value. */
  itkSetMacro( PaddingValue, FixedImagePixelType );

  /** Returns the padding value. */
  itkGetConstReferenceMacro( PaddingValue, FixedImagePixelType );


  /** Return the joint histogram. This is updated during every call to the 
   *  GetValue() method. The histogram can for instance be used by 
   *  itk::HistogramToImageFilter to plot the joint histogram. */
  itkGetConstReferenceMacro( Histogram, HistogramPointer );
  
  /** Set whether the padding value should be used to determine which pixels
      should be ignored when calculating the similarity measure. Those pixels
      in the fixed image which have the padding value will be ignored. */
  itkSetMacro( UsePaddingValue, bool );
  itkGetConstMacro( UsePaddingValue, bool );

  //itkSetMacro( UseAllPixels, bool );
  //itkGetConstMacro( UseAllPixels, bool );

  itkSetMacro( UseCoincidenceWeighting, bool );
  itkGetConstMacro( UseCoincidenceWeighting, bool );
  itkSetMacro( WeightingFactor, double );
  itkGetConstMacro( WeightingFactor, double );
  itkSetMacro( UseAdaptativeBining, bool );
  itkGetConstMacro( UseAdaptativeBining, bool );
  

  /** Sets the step length used to calculate the derivative. */
  itkSetMacro( DerivativeStepLength, double );

  /** Returns the step length used to calculate the derivative. */
  itkGetConstMacro( DerivativeStepLength, double );

  /** The scales type. */
  typedef Array<double> ScalesType;

  /** Sets the derivative step length scales. */
  itkSetMacro( DerivativeStepLengthScales, ScalesType );

  /** Returns the derivate step length scales. */
  itkGetConstReferenceMacro(DerivativeStepLengthScales, ScalesType);

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType& parameters) const;

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters,
                     DerivativeType & derivative) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters,
                             MeasureType& Value,
                             DerivativeType& Derivative) const;

  /** Set the lower bounds of the intensities to be considered for computing
    * the histogram. This option allows to focus the computation of the Metric in
    * a particular range of intensities that correspond to features of interest. */
  void SetLowerBound( const MeasurementVectorType & bound );

  /** Set the upper bounds of the intensities to be considered for computing
    * the histogram. This option allows to focus the computation of the Metric in
    * a particular range of intensities that correspond to features of interest.  */
  void SetUpperBound( const MeasurementVectorType & bound );


  void SetBackgroundThresholds( const MeasurementVectorType & thres );


  itkGetConstReferenceMacro( BinThresholds, HistogramType::IndexType  );
  itkGetConstReferenceMacro( Thresholds, MeasurementVectorType  );
  
  itkSetMacro( NumberOfSpatialSamples, unsigned int );
  itkGetConstReferenceMacro( NumberOfSpatialSamples, unsigned int );
protected:
  /** Constructor is protected to ensure that \c New() function is used to
      create instances. */
  HistogramImageToImageMetric();
  virtual ~HistogramImageToImageMetric() {};

  /** The histogram size. */
  HistogramSizeType m_HistogramSize;
  /** The lower bound for samples in the histogram. */
  mutable MeasurementVectorType m_LowerBound;
  /** The upper bound for samples in the histogram. */
  mutable MeasurementVectorType m_UpperBound;

  /** The maxs container for samples in the histogram */
  BinMaxContainerType m_Max;
  BinMinContainerType m_Min;

  /** The increase in the upper bound. */
  double m_UpperBoundIncreaseFactor;

  /** Boolean flag to indicate whether the user supplied lower bounds or
    * whether they should be computed from the min of image intensities */
  bool m_LowerBoundSetByUser;

  /** Boolean flag to indicate whether the user supplied upper bounds or
    * whether they should be computed from the max of image intensities */
  bool m_UpperBoundSetByUser;

  MeasurementVectorType m_Thresholds;
  mutable typename HistogramType::IndexType m_BinThresholds;
  bool m_ThresholdsSetByUser;
  mutable double m_WeightingFactor;
  bool m_UseCoincidenceWeighting;
  bool m_UseAdaptativeBining;

  /** Computes the joint histogram from the transformation parameters
      passed to the function. */
  void ComputeHistogram(const TransformParametersType & parameters,
                        HistogramType& histogram) const;
  /** Computes the joint histogram from the transformation parameters
      passed to the function. */
  void ComputeHistogram(const TransformParametersType & parameters,
                        unsigned int parameter,
                        double step,
                        HistogramType& histogram) const;
  /** Copies a histogram. */
  void CopyHistogram(HistogramType& target, HistogramType& source) const;

  /** Evaluates the similarity measure using the given histogram. All
      subclasses must reimplement this method. */
  virtual MeasureType EvaluateMeasure(HistogramType& histogram) const = 0;
  
  /** PrintSelf funtion */
  void PrintSelf(std::ostream& os, Indent indent) const;

  void FillSubsampledHistogram( HistogramType& histogram ) const;

private:
  HistogramImageToImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The padding value. */
  FixedImagePixelType m_PaddingValue;

  /** True if those pixels in the fixed image with the same value as the
      padding value should be ignored when calculating the similarity
      measure. */
  bool m_UsePaddingValue;
  //bool m_UseAllPixels;

  /** The step length used to calculate the derivative. */
  double m_DerivativeStepLength;

  unsigned int m_NumberOfSpatialSamples;

  /** The derivative step length scales. */
  ScalesType m_DerivativeStepLengthScales;

  /** Pointer to the joint histogram. This is updated during every call to 
   * GetValue() */
  HistogramPointer  m_Histogram;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHistogramImageToImageMetric.txx"
#endif

#endif // __itkHistogramImageToImageMetric_h
