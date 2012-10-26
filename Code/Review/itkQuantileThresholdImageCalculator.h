/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuantileThresholdImageCalculator.h,v $
  Language:  C++
  Date:      $Date: 2009-04-23 03:53:36 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuantileThresholdImageCalculator_h
#define __itkQuantileThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

#include "../Review/itkHistogram.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkSampleToHistogramFilter.h"

namespace itk
{

/** \class QuantileThresholdImageCalculator
 * \brief Computes the Otsu's threshold for an image.
 * 
 * This calculator computes the Otsu's threshold which separates an image
 * into foreground and background components. The method relies on a
 * histogram of image intensities. The basic idea is to maximize the 
 * between-class variance.
 *
 * This class is templated over the input image type.
 *
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT QuantileThresholdImageCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef QuantileThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(QuantileThresholdImageCalculator, Object);

  /** Type definition for the input image. */
  typedef TInputImage  ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer  ImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;
  
  /** Type definition for the input image region type. */
  typedef typename TInputImage::RegionType RegionType;
  
#ifdef ITK_USE_REVIEW_STATISTICS
  typedef Statistics::Histogram<double>                  HistogramType;
#else
  typedef Statistics::Histogram<double, 1>               HistogramType;
#endif

  typedef typename HistogramType::MeasurementVectorType  MeasurementVectorType;
  typedef typename HistogramType::SizeType               HistogramSizeType;
  typedef typename HistogramType::Pointer                HistogramPointer;
  typedef typename HistogramType::BinMinContainerType    BinMinContainerType;
  typedef typename HistogramType::BinMaxContainerType    BinMaxContainerType;

  typedef itk::Statistics::ImageToListSampleAdaptor<TInputImage>     ImageAdaptorType;
  typedef itk::Statistics::SampleToHistogramFilter< ImageAdaptorType, HistogramType >    HistogramFilter;

  /** Set the input image. */
  itkSetConstObjectMacro(Image,ImageType);

  /** Compute the Otsu's threshold for the input image. */
  void Compute(void);

  /** Return the Otsu's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1,
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  itkSetClampMacro( Percentage, double, 0.0, 1.0 );
  itkGetConstMacro( Percentage, double );

protected:
  QuantileThresholdImageCalculator();
  virtual ~QuantileThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  QuantileThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  ImageConstPointer    m_Image;
  double               m_Percentage;
  unsigned long        m_NumberOfHistogramBins;
  HistogramType::ConstPointer     m_Histogram;
  bool                 m_HistogramComputed;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuantileThresholdImageCalculator.txx"
#endif

#endif
