/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSampleToAnisotropicHistogramFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-05-02 05:43:58 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSampleToAnisotropicHistogramFilter_h

#ifndef _WIN32
	#warning "Using SampleToAnisotropicHistogramFilter"
#endif

#define __itkSampleToAnisotropicHistogramFilter_h

#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkMeasurementVectorTraits.h"
#include "itkSimpleDataObjectDecorator.h"
#include "itkSampleToHistogramFilter.h"


namespace itk { 
namespace Statistics {

/** \class SampleToAnisotropicHistogramFilter 
 *  \brief Computes the Histogram corresponding to a Sample.
 *
 * This filter produces as output the histogram corresponding to
 * the values of a Sample.
 *
 * \sa Sample, Histogram
 *
 */

template < class TSample, class THistogram >
class ITK_EXPORT SampleToAnisotropicHistogramFilter : public ProcessObject
{
public:
  /** Standard class typedefs */
  typedef SampleToAnisotropicHistogramFilter       Self;  
  typedef ProcessObject                 Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer<const Self >     ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(SampleToAnisotropicHistogramFilter, ProcessObject);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** MeasurementVector typedef support */ 
  typedef TSample                                        SampleType;
  typedef THistogram                                     HistogramType;
  typedef typename SampleType::MeasurementVectorType     MeasurementVectorType;
  typedef typename MeasurementVectorType::ValueType      MeasurementType;
  typedef typename HistogramType::SizeType               HistogramSizeType;
  typedef typename HistogramType::MeasurementType        HistogramMeasurementType;
  typedef typename HistogramType::MeasurementVectorType  HistogramMeasurementVectorType;


  /** Type for the data object output */
  itkSuperclassTraitMacro( DataObjectPointer )

  /** Set/Get the input sample */
  virtual void SetInput( const SampleType * sample );
  virtual const SampleType * GetInput() const;
  
  /** Get the output Histogram */
  const HistogramType  * GetOutput() const;

  /** Type of DataObjects to use for Size inputs */
  typedef SimpleDataObjectDecorator<
    HistogramSizeType > InputHistogramSizeObjectType;

  /** Type of DataObjects to use for Marginal Scale inputs */
  typedef SimpleDataObjectDecorator<
    HistogramMeasurementType > InputHistogramMeasurementObjectType;

  /** Type of DataObjects to use for Minimum and Maximums values of the
   * histogram bins. */
  typedef SimpleDataObjectDecorator<
    HistogramMeasurementVectorType > InputHistogramMeasurementVectorObjectType;

  /** Type of DataObjects to use for AutoMinimumMaximum input */
  typedef SimpleDataObjectDecorator< bool > InputBooleanObjectType;

  typedef SampleToHistogramFilter
	       < SampleType, HistogramType > PDFHistogramFilter;

  /** Methods for setting and getting the histogram size.  The histogram size
   * is encapsulated inside a decorator class. For this reason, it is possible
   * to set and get the decorator class, but it is only possible to set the
   * histogram size by value. This macro declares the methods
   * SetHistogramSize(), SetHistogramSizeInput(), GetHistogramSizeInput().
   */
  itkSetDecoratedInputMacro( HistogramSize, HistogramSizeType, 1 );

  /** Methods for setting and getting the Marginal scale value.  The marginal
   * scale is used when the type of the measurement vector componets are of
   * integer type. */
  itkSetDecoratedInputMacro( MarginalScale, HistogramMeasurementType, 2 );

  /** Methods for setting and getting the Minimum and Maximum values of the
   * histogram bins. */
  itkSetDecoratedInputMacro( HistogramBinMinimum, HistogramMeasurementVectorType, 3 );
  itkSetDecoratedInputMacro( HistogramBinMaximum, HistogramMeasurementVectorType, 4 );

  /** Methods for setting and getting the boolean flag that defines whether the
   * minimum and maximum of the histogram are going to be computed
   * automatically from the values of the sample */
  itkSetDecoratedInputMacro( AutoMinimumMaximum, bool, 5 );

  /** Method for setting the density factor **/
  itkSetDecoratedInputMacro( PDFDensityFactor, unsigned int, 6 );

  itkSetDecoratedInputMacro( MaxEstimatorIterations, unsigned int, 7 );


  /** Method that facilitates the use of this filter in the internal
   * pipeline of another filter. */
  virtual void GraftOutput(DataObject *output);

  double GetNoisePower ( unsigned int dim );
  double GetSignalPower( unsigned int dim );
  double GetErrorSum   ( unsigned int dim );

protected:
  SampleToAnisotropicHistogramFilter();
  virtual ~SampleToAnisotropicHistogramFilter();

  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Make a DataObject of the correct type to used as the specified
   * output. This method
   * is automatically called when DataObject::DisconnectPipeline() is
   * called.  
   * \sa ProcessObject
   */
  virtual DataObjectPointer MakeOutput(unsigned int idx);
  
  // Where the histogram is actually computed
  virtual void GenerateData();

  virtual void Sample();

  virtual void GeneratePDF(HistogramMeasurementVectorType lower, HistogramMeasurementVectorType upper);

  void ComputeBinWidths( unsigned int dim );

private:
  SampleToAnisotropicHistogramFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename PDFHistogramFilter::Pointer pdfGenerator;
  const HistogramType *  m_pdf;

}; // end of class

} // end of namespace Statistics 
} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSampleToAnisotropicHistogramFilter.txx"
#endif
  
#endif
