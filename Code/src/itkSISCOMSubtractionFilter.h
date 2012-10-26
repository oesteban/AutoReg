/* --------------------------------------------------------------------------------------
 * File:    itkSISCOMSubtractionFilter.h
 * Date:    01/08/2011
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


#ifndef ITKSISCOMSUBTRACTIONFILTER_H_
#define ITKSISCOMSUBTRACTIONFILTER_H_

#include "itkImageToImageFilter.h"
//#include "itkFixedArray.h"
#include "itkImage.h"
//#include "itkSubtractImageFilter.h"
#include "itkGIBUBSubtractionFilter.h"
#include "itkMeanImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkSpatialObject.h"

namespace itk
{
/**
 * \class SISCOMSubtractionFilter
 * \brief Class for ictal and interictal subtraction on
 * SISCOM methodology
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 *
 * \ingroup ImageEnhancement
 * \ingroup ImageFeatureExtraction
 */

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType=double >
class ITK_EXPORT SISCOMSubtractionFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef SISCOMSubtractionFilter                         Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SISCOMSubtractionFilter, ImageToImageFilter);

  /** Image type information. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType         OutputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputInternalPixelType;
  typedef typename TOutputImage::Pointer           OutputImagePointer;
  typedef typename TOutputImage::ConstPointer      OutputImageConstPointer;

  typedef typename TInputImage::PixelType          InputPixelType;
  typedef typename TInputImage::InternalPixelType  InputInternalPixelType;
  typedef typename TInputImage::Pointer            InputImagePointer;
  typedef typename TInputImage::ConstPointer       InputImageConstPointer;

  typedef itk::MultiplyByConstantImageFilter<InputImageType, double, InputImageType> ApplyFactorFilter;
  typedef typename itk::ImageRegionConstIterator< InputImageType > InputImageConstIterator;
  typedef typename itk::ImageRegionIterator< OutputImageType >     OutputImageIterator;

  itkStaticConstMacro(InputImageDimension, unsigned int, TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int, TOutputImage::ImageDimension);

  typedef itk::MinimumMaximumImageCalculator <InputImageType>	                           CalculatorFilter;
  typedef itk::GIBUBSubtractionFilter < InputImageType, InputImageType , OutputImageType > SubtractFilter;
  typedef itk::MeanImageFilter< InputImageType, InputImageType >                           MeanFilter;
  typedef itk::ResampleImageFilter< InputImageType, InputImageType>                        ResampleFilter;

  //  typedef itk::SubtractImageFilter< InputImageType, InputImageType, OutputImageType> SubtractFilter;

  /** Typedef of double containers */
  //typedef FixedArray<double, itkGetStaticConstMacro(ImageDimension)> ArrayType;

  itkSetMacro(UseNormalization,bool);
  itkGetConstMacro(UseNormalization, const bool);

  /** Transform typedef. */
  typedef MatrixOffsetTransformBase<TInterpolatorPrecisionType,
    itkGetStaticConstMacro(InputImageDimension),
    itkGetStaticConstMacro(InputImageDimension)> TransformType;
  typedef typename TransformType::Pointer      TransformPointerType;
  typedef typename TransformType::ConstPointer TransformConstPointerType;

  /** Interpolator typedef. */
  typedef InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType> InterpolatorType;
  typedef typename InterpolatorType::Pointer  InterpolatorPointerType;


  /**  Type for the mask of the moving image. Only pixels that are "inside"
       this mask will be considered for the computation of the subtraction */
  typedef SpatialObject< itkGetStaticConstMacro(InputImageDimension) >   InputImageMaskType;
  typedef typename InputImageMaskType::Pointer                           InputImageMaskPointer;
  typedef typename InputImageMaskType::ConstPointer                      InputImageMaskConstPointer;

  /** base type for images of the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(InputImageDimension)> ImageBaseType;

  /** Set the coordinate transformation.
   * Set the coordinate transform to use for resampling the interictal image to
   * ictal image space.
   * By default the filter uses an Identity transform. You must provide a different
   * transform here, before attempting to run the filter, if you do not want to
   * use the default Identity transform. */
  itkSetConstObjectMacro( RegTransform, TransformType );

  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( RegTransform, TransformType );

  /** Set the coordinate transformation.
   * Set the coordinate transform to use for resampling the interictal image to
   * RM reference image space.
   * By default the filter uses an Identity transform. You must provide a different
   * transform here, before attempting to run the filter, if you do not want to
   * use the default Identity transform. */
  itkSetConstObjectMacro( InterictalTransform, TransformType );
  /** Get a pointer to the coordinate transform. */
  itkGetConstObjectMacro( InterictalTransform, TransformType );

  /** Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>. Some
   * other options are itk::NearestNeighborInterpolateImageFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and itk::BSplineInterpolateImageFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro( DefaultPixelValue, InputPixelType );

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetConstReferenceMacro( DefaultPixelValue, InputPixelType );


  itkGetConstReferenceMacro( NormalizationFactor, double );

  /** Connect one of the operands for pixel-wise addition */
  void SetInterictalImage( const TInputImage * image);

  /** Connect one of the operands for pixel-wise addition */
  void SetIctalImage( const TInputImage * image);

  /** Set/Get the interictal image mask. */
  itkSetObjectMacro( InterictalImageMask, InputImageMaskType );

#ifdef ITK_LEGACY_REMOVE
  itkSetConstObjectMacro( InterictalImageMask, InputImageMaskType );
#else
  virtual void SetInterictalImageMask( const InputImageMaskType* mask )
    { this->SetInterictalImageMask(const_cast<InputImageMaskType*>(mask)); }
#endif
  itkGetConstObjectMacro( InterictalImageMask, InputImageMaskType );

  /** Set/Get the ictal image mask. */
  itkSetObjectMacro( IctalImageMask, InputImageMaskType );
#ifdef ITK_LEGACY_REMOVE
  itkSetConstObjectMacro( IctalImageMask, InputImageMaskType );
#else
  virtual void SetIctalImageMask( const InputImageMaskType* mask )
    { this->SetIctalImageMask(const_cast<InputImageMaskType*>(mask)); }
#endif
  itkGetConstObjectMacro( IctalImageMask, InputImageMaskType );


  /** Get the image output of this process object.  */
  InputImageConstPointer GetInterictalNormalizedOutput() { return  m_InputSource; }

  itkGetConstObjectMacro( IctalRefSpace, InputImageType );
  itkGetConstObjectMacro( InterictalRefSpace, InputImageType );
  itkGetConstObjectMacro( SubtractionRefSpace, OutputImageType );

  TransformPointerType GetIctalTransform();
  void ResampleToReferenceSpace( const InputImageType* i_ref );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<OutputPixelType>));

  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
                            itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif

protected:
  SISCOMSubtractionFilter();
  virtual ~SISCOMSubtractionFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Standard pipeline method. While this class does not implement a
   * ThreadedGenerateData(), its GenerateData() delegates all
   * calculations to an NeighborhoodOperatorImageFilter.  Since the
   * NeighborhoodOperatorImageFilter is multithreaded, this filter is
   * multithreaded by default. */
  void GenerateData();


private:
  SISCOMSubtractionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void ComputeNormalizationFactor();
  void ComposeTransforms();

  /** Flag to indicate whether to normalize input 1 */
  bool m_UseNormalization;

  double m_NormalizationFactor;

  InputPixelType               m_DefaultPixelValue;
  InputImageConstPointer       m_InputSource;
  TransformConstPointerType    m_RegTransform;         // Coordinate transform to use
  TransformConstPointerType    m_InterictalTransform;         // Coordinate transform to use
  TransformPointerType         m_IctalTransform;         // Coordinate transform to use
  InterpolatorPointerType      m_Interpolator;

  InputImageConstPointer       m_IctalRefSpace;
  InputImageConstPointer       m_InterictalRefSpace;
  OutputImageConstPointer      m_SubtractionRefSpace;

  typename SubtractFilter::Pointer m_Subtract;
  typename ResampleFilter::Pointer m_Resample;

#ifdef ITK_LEGACY_REMOVE
  InputImageMaskConstPointer m_InterictalImageMask;
  InputImageMaskConstPointer m_IctalImageMask;
#else
  mutable InputImageMaskPointer m_InterictalImageMask;
  mutable InputImageMaskPointer m_IctalImageMask;
#endif

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSISCOMSubtractionFilter.txx"
#endif


#endif /* ITKSISCOMSUBTRACTIONFILTER_H_ */
