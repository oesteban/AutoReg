/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHistogram.h,v $
  Language:  C++
  Date:      $Date: 2009-05-22 12:54:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkHistogram_h

#ifndef _WIN32
	#warning "Using ReviewHistogram"
#endif

#define __itkHistogram_h

#include <vector>

#include "itkArray.h"
#include "itkSample.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkSparseFrequencyContainer2.h"

namespace itk {
namespace Statistics {

/**
 * Due to a bug in MSVC, an enum value cannot be accessed out of a template
 * parameter until the template class opens.  In order for templated classes
 * to access the dimension of an image template parameter in defining their
 * own dimension, this class is needed as a work-around.
 */
template <typename THistogram>
struct GetHistogramDimension
{
  itkStaticConstMacro(HistogramDimension, unsigned int, THistogram::MeasurementVectorSize);
};

/** \class Histogram
 *  \brief This class stores measurement vectors in the context of n-dimensional histogram.
 *
 * Histogram represents an ND histogram.  Histogram bins can be
 * regularly or irregularly spaced. The storage for the histogram is
 * managed via the FrequencyContainer specified by the template
 * argument.  The default frequency container is a
 * DenseFrequencyContainer. A SparseFrequencyContainer can be used as
 * an alternative.
 *
 * Frequencies of a bin (SetFrequency(), IncreaseFrequency()) can be
 * specified by measurement, index, or instance identifier.
 *
 * Measurements can be queried by bin index or instance
 * identifier. In this case, the measurement returned in the centroid
 * of the histogram bin.
 *
 * The Initialize() method is used to specified the number of bins for
 * each dimension of the histogram. An overloaded version also allows
 * for regularly spaced bins to defined.  To define irregularly sized
 * bins, use the SetBinMin()/SetBinMax() methods.
 *
 * If you do not know the length of the measurement vector at compile time, you
 * should use the VariableDimensionHistogram class, instead of the Histogram
 * class.
 *
 * If you know the length of the measuremen vector at compile time, you can
 * conveniently be obtained from MeasurementVectorTraits. For instance,
 * instantiate a histogram as below:
 *
 * \code
 * typedef Histogram< THistogramMeasurement, typename TFrequencyContainer > HistogramType;
 * \endcode
 *
 * \sa Sample, DenseFrequencyContainer, SparseFrequencyContainer, VariableDimensionHistogram
 */

template < class TMeasurement = float,
           class TFrequencyContainer = DenseFrequencyContainer2 >
class ITK_EXPORT Histogram
  : public Sample < Array< TMeasurement > >
{
public:

  // This type serves as the indirect definition of MeasurementVectorType
  typedef Array< TMeasurement >        ArrayType;

  /** Standard typedefs */
  typedef Histogram                    Self;
  typedef Sample< ArrayType  >         Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(Histogram, Sample);

  /** standard New() method support */
  itkNewMacro(Self);

  /** type of an element of a measurement vector */
  typedef TMeasurement MeasurementType;

  /** Common sample class typedefs */
  itkSuperclassTraitMacro(MeasurementVectorType)
  itkSuperclassTraitMacro(InstanceIdentifier)
  itkSuperclassTraitMacro(MeasurementVectorSizeType)
    typedef MeasurementVectorType ValueType;

  /** frequency container typedef */
  typedef TFrequencyContainer                        FrequencyContainerType;
  typedef typename FrequencyContainerType::Pointer   FrequencyContainerPointer;

  /** Frequency and TotalFrequency value type from superclass */
  typedef typename FrequencyContainerType::AbsoluteFrequencyType      AbsoluteFrequencyType;
  typedef typename FrequencyContainerType::TotalAbsoluteFrequencyType TotalAbsoluteFrequencyType;
  typedef typename FrequencyContainerType::RelativeFrequencyType      RelativeFrequencyType;
  typedef typename FrequencyContainerType::TotalRelativeFrequencyType TotalRelativeFrequencyType;
  
  /** Index typedef support. An index is used to access pixel values. */
  typedef Array< long >                     IndexType;
  typedef typename IndexType::ValueType     IndexValueType;

  /** size array type */
  typedef Array< unsigned long >            SizeType;
  typedef typename SizeType::ValueType      SizeValueType;

  /** bin min max value storage types */
  typedef std::vector< MeasurementType >    BinMinVectorType;
  typedef std::vector< MeasurementType >    BinMaxVectorType;
  typedef std::vector< BinMinVectorType >   BinMinContainerType;
  typedef std::vector< BinMaxVectorType >   BinMaxContainerType;


  /** Initialize the histogram, generating the offset table and
   * preparing the frequency container. Subclasses should call this
   * method in their Initialize() method. */
  void Initialize(const SizeType &size);


  /** Initialize the histogram using equal size bins. To assign bin's
   * min and max values along each dimension use SetBinMin() and
   * SetBinMax() functions. */
  void Initialize(const SizeType &size, MeasurementVectorType& lowerBound,
                  MeasurementVectorType& upperBound);

  /** Initialize the histogram using explicitly passed bin's
   * min and max values along each dimension use SetBinMin() and
   * SetBinMax() functions. */
  void Initialize(const BinMinContainerType& mins, const BinMaxContainerType& maxs);

  /** Initialize the values of the histogram bins to zero */
  void SetToZero();

  /** Get the index of a measurement value from the histogram.
   * \deprecated Use GetIndex(const MeasurementVectorType &
   * measurement, IndexType & index ) const instead. */
  const IndexType & GetIndex(const MeasurementVectorType& measurement) const;

  /** Get the index of histogram corresponding to the specified
   *  measurement value. Returns true if index is valid and false if
   *  the measurement is outside the histogram */
  bool GetIndex(const MeasurementVectorType & measurement,
                IndexType & index ) const;

  /** Get the index that is uniquely labelled by an instance identifier
   * The corresponding id is the offset of the index
   * This method uses ImageBase::ComputeIndex() method */
  const IndexType & GetIndex( InstanceIdentifier id) const;

  /** Is set to false if the bins at edges of the histogram extend to
   *   +/- infinity. */
  itkGetConstMacro(ClipBinsAtEnds, bool);

  /** Set to false to have the bins at edges of the histogram extend to
   *   +/- infinity. */
  itkSetMacro(ClipBinsAtEnds, bool);

  /** Returns true if the given index is out of bound meaning one of index
   * is not between [0, last index] */
  bool IsIndexOutOfBounds(const IndexType &index) const;

  /** Get the instance identifier of the bin that is indexed by the
   * index. The corresponding instance identifier is the offset of the index
   * This method uses ImageBase::ComputeIndex() method */
  InstanceIdentifier GetInstanceIdentifier(const IndexType &index) const;

  /** Returns the number of instances (bins or cells) in this container */
  InstanceIdentifier Size() const;

  /** Get the size (N-dimensional) of the histogram  */
  const SizeType & GetSize() const;

  /** Get the size of histogram along a specified dimension */
  SizeValueType GetSize(unsigned int dimension) const;

  /** Get the minimum value of nth bin of dimension d */
  const MeasurementType& GetBinMin(unsigned int dimension,
                             InstanceIdentifier nbin) const;

  /** Get the maximum value of nth bin of dimension d */
  const MeasurementType& GetBinMax(unsigned int dimension,
                             InstanceIdentifier nbin) const;

  /** Set the minimum value of nth bin of dimension d */
  void SetBinMin(unsigned int dimension, InstanceIdentifier nbin,
                 MeasurementType min);

  /** Set the maximum value of nth bin of dimension d */
  void SetBinMax(unsigned int dimension,
                 InstanceIdentifier nbin, MeasurementType max);

  /** Get the minimum of the bin along dimension d corresponding to a
   * particular measurement. */
  const MeasurementType& GetBinMinFromValue(unsigned int dimension,
                                      float value ) const;

  /** Get the maximum of the bin along dimension d corresponding to a
   * particular measurement. */
  const MeasurementType& GetBinMaxFromValue(unsigned int dimension,
                                      float value ) const;

  /** Get the vector of bin minimums along a dimension  */
  const BinMinVectorType& GetDimensionMins(unsigned int dimension) const;

  /** Get the vector of maximums along a dimension  */
  const BinMaxVectorType& GetDimensionMaxs(unsigned int dimension) const;

  /** Get the minimums of the bins  */
  const BinMinContainerType& GetMins() const;

  /** Method the maximums of the bins  */
  const BinMaxContainerType& GetMaxs() const;

  /** Get the minimums of the bin corresponding to a particular index */
  const MeasurementVectorType & GetHistogramMinFromIndex(const IndexType &index) const;

  /** Get the maximums of the bin corresponding to a particular index  */
  const MeasurementVectorType& GetHistogramMaxFromIndex(const IndexType &index) const;

  /** Get the frequency of an instance indentifier */
  AbsoluteFrequencyType GetFrequency( InstanceIdentifier id ) const;

  /** Get the frequency of an index */
  AbsoluteFrequencyType GetFrequency(const IndexType &index) const;

  /** Set all the bins in the histogram to a specified frequency */
  void SetFrequency( AbsoluteFrequencyType value );

  /** Set the frequency of an instance identifier.  Returns false if the bin is
   * out of bounds. */
  bool SetFrequency( InstanceIdentifier id, AbsoluteFrequencyType value);

  /** Set the frequency of an index. Returns false if the bin is
   * out of bounds. */
  bool SetFrequency(const IndexType &index,
                    AbsoluteFrequencyType value);

  /** Set the frequency of a measurement. Returns false if the bin is
   * out of bounds. */
  bool SetFrequency(const MeasurementVectorType &measurement,
                    AbsoluteFrequencyType value);


  /** Increase the frequency of an instance identifier.
   * Frequency is increased by the specified value. Returns false if
   * the bin is out of bounds. */
  bool IncreaseFrequency(InstanceIdentifier id, AbsoluteFrequencyType value);

  /** Increase the frequency of an index.  Frequency is
   * increased by the specified value. Returns false if the bin is out
   * of bounds. */
  bool IncreaseFrequency(const IndexType &index, AbsoluteFrequencyType value);

  /** Increase the frequency of a measurement.  Frequency is
   * increased by the specified value. Returns false if the
   * measurement is outside the bounds of the histogram. */
  bool IncreaseFrequency(const MeasurementVectorType &measurement,
                         AbsoluteFrequencyType value);

  /** Get the measurement of an instance identifier. This is the
   * centroid of the bin.
   */
  const MeasurementVectorType & GetMeasurementVector( InstanceIdentifier id) const;

  /** Get the measurement of an index. This is the centroid of the bin. */
  const MeasurementVectorType & GetMeasurementVector(const IndexType &index) const;

  /** Get the measurement a bin along a specified dimension.  This is
   * the midpoint of the bin along that dimension. */
  MeasurementType GetMeasurement(InstanceIdentifier n,
                                 unsigned int dimension) const;

  /** Get the total frequency in the histogram */
  TotalAbsoluteFrequencyType GetTotalFrequency() const;

  /** Get the frequency of a dimension's nth element. */
  AbsoluteFrequencyType GetFrequency(InstanceIdentifier n,
                             unsigned int dimension) const;

  /** Get the pth percentile value for a dimension.
   *
   * Let assume n = the index of the bin where the p-th percentile value is,
   * min = min value of the dimension of the bin,
   * max = max value of the dimension of the bin,
   * interval = max - min ,
   * pp = cumlated proportion until n-1 bin;
   * and pb = frequency of the bin / total frequency of the dimension.
   *
   * If p is less than 0.5,
   * the percentile value =
   * min + ((p - pp ) / pb) * interval
   * If p is greater than or equal to 0.5
   * the percentile value =
   * max - ((pp - p) / pb) * interval  */
  double Quantile(unsigned int dimension, double p) const;

  /** Method to graft another histogram's output */
  virtual void Graft( const DataObject * );

protected:
  void PrintSelf(std::ostream& os, Indent indent) const;
  void PrintMaxs(std::ostream& os, Indent indent) const;
  void PrintMins(std::ostream& os, Indent indent) const;

public:

  /** \class Histogram::ConstIterator class that walks through the elements of
   * the histogram */
  class ConstIterator
    {
    public:

    friend class Histogram;

    ConstIterator(const Self * histogram)
      {
      m_Id = 0;
      m_Histogram = histogram;
      }

    ConstIterator(const ConstIterator & it)
      {
      m_Id        = it.m_Id;
      m_Histogram = it.m_Histogram;
      }

    ConstIterator& operator=(const ConstIterator& it)
      {
      m_Id  = it.m_Id;
      m_Histogram = it.m_Histogram;
      return *this;
      }

    AbsoluteFrequencyType GetFrequency() const
      {
      return  m_Histogram->GetFrequency(m_Id);
      }

    InstanceIdentifier GetInstanceIdentifier() const
      {
      return m_Id;
      }

    const MeasurementVectorType & GetMeasurementVector() const
      {
      return m_Histogram->GetMeasurementVector(m_Id);
      }

    ConstIterator& operator++()
      {
      ++m_Id;
      return *this;
      }

    bool operator!=(const ConstIterator& it)
      {
      return (m_Id != it.m_Id);
      }

    bool operator==(const ConstIterator& it)
      {
      return (m_Id == it.m_Id);
      }

#if !(defined(_MSC_VER) && (_MSC_VER <= 1200))
  protected:
#endif
    // This method is purposely not implemented
    ConstIterator();

    ConstIterator(InstanceIdentifier id, const Self * histogram)
      : m_Id(id), m_Histogram(histogram)
      {}

    // ConstIterator pointing DenseFrequencyContainer
    InstanceIdentifier m_Id;

    // Pointer of DenseFrequencyContainer
    const Self* m_Histogram;
    }; // end of iterator class

  /** \class Histogram::Iterator */
  class Iterator : public ConstIterator
  {
  public:

    Iterator(Self * histogram):ConstIterator( histogram )
      {
      }

    Iterator(InstanceIdentifier id, Self * histogram)
        : ConstIterator( id, histogram )
      {}

    Iterator(const Iterator & it):ConstIterator(it)
      {
      }

    Iterator & operator=(const Iterator& it)
      {
      this->ConstIterator::operator=( it );
      return *this;
      }

    bool SetFrequency(const AbsoluteFrequencyType value)
      {
      Self * histogram = const_cast< Self * >( this->m_Histogram );
      return histogram->SetFrequency( this->m_Id, value );
      }

#if !(defined(_MSC_VER) && (_MSC_VER <= 1200))
  protected:
#endif
    // To ensure const-correctness these method must not be in the public API.
    // The are purposly not implemented, since they should never be called.
    Iterator();
    Iterator(const Self * histogram);
    Iterator(InstanceIdentifier id, const Self * histogram);
    Iterator(const ConstIterator & it);
    ConstIterator& operator=(const ConstIterator& it);

  private:
    }; // end of iterator class


  Iterator  Begin()
    {
    Iterator iter(0, this);
    return iter;
    }

  Iterator  End()
    {
    return Iterator(m_OffsetTable[this->GetMeasurementVectorSize()], this);
    }

  ConstIterator  Begin() const
    {
    ConstIterator iter(0, this);
    return iter;
    }

  ConstIterator End() const
    {
    return ConstIterator(m_OffsetTable[this->GetMeasurementVectorSize()], this);
    }

protected:
  Histogram();
  virtual ~Histogram() {}

  // The number of bins for each dimension
  SizeType m_Size;

private:
  Histogram(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typedef std::vector< InstanceIdentifier >   OffsetTableType;
  OffsetTableType                             m_OffsetTable;
  FrequencyContainerPointer                   m_FrequencyContainer;
  unsigned int                                m_NumberOfInstances;

  // This method is provided here just to avoid a "hidden" warning 
  // related to the virtual method available in DataObject.
  virtual void Initialize() {};

  // lower bound of each bin
  std::vector< std::vector<MeasurementType> > m_Min;

  // upper bound of each bin
  std::vector< std::vector<MeasurementType> > m_Max;

  mutable MeasurementVectorType   m_TempMeasurementVector;
  mutable IndexType               m_TempIndex;

  bool                            m_ClipBinsAtEnds;
};

} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHistogram.txx"
#endif

#endif
