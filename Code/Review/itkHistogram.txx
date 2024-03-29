/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHistogram.txx,v $
  Language:  C++
  Date:      $Date: 2009-05-04 15:47:43 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkHistogram_txx
#define __itkHistogram_txx

#include "itkHistogram.h"
#include "itkNumericTraits.h"

namespace itk { 
namespace Statistics {

template< class TMeasurement, class TFrequencyContainer >
Histogram<TMeasurement, TFrequencyContainer >
::Histogram()
{
  this->m_ClipBinsAtEnds = true;
  this->m_FrequencyContainer = FrequencyContainerType::New();
  this->m_OffsetTable = OffsetTableType( this->GetMeasurementVectorSize() + 1 );
  for( unsigned int i = 0; i < this->GetMeasurementVectorSize() + 1; i++ )
    {
    this->m_OffsetTable[i] = itk::NumericTraits< InstanceIdentifier >::Zero;
    }
}

template< class TMeasurement, class TFrequencyContainer >
typename 
Histogram<TMeasurement, TFrequencyContainer >::InstanceIdentifier
Histogram<TMeasurement, TFrequencyContainer >
::Size() const
{
  if( this->GetMeasurementVectorSize() == 0 )
    {
    return itk::NumericTraits< InstanceIdentifier >::Zero;
    }
  InstanceIdentifier size = 1;
  for (unsigned int i = 0; i < this->GetMeasurementVectorSize(); i++)
    {
    size *= m_Size[i];
    }
  return size;
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::SizeType &
Histogram<TMeasurement, TFrequencyContainer >
::GetSize() const
{
  return m_Size;
}

template< class TMeasurement, class TFrequencyContainer >
typename 
Histogram<TMeasurement, TFrequencyContainer >::SizeValueType
Histogram<TMeasurement, TFrequencyContainer >
::GetSize(unsigned int dimension) const
{
  return m_Size[dimension];
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::MeasurementType &
Histogram<TMeasurement, TFrequencyContainer >
::GetBinMin(unsigned int dimension, InstanceIdentifier nbin) const
{
  return m_Min[dimension][nbin];
}

template< class TMeasurement, 
          class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::MeasurementType &
Histogram<TMeasurement, TFrequencyContainer >
::GetBinMax(unsigned int dimension, InstanceIdentifier nbin) const
{
  return m_Max[dimension][nbin];
}

template< class TMeasurement, class TFrequencyContainer >
void
Histogram<TMeasurement, TFrequencyContainer >
::SetBinMin( unsigned int dimension, InstanceIdentifier nbin,
                 MeasurementType min)
{
  m_Min[dimension][nbin] = min;
}

template< class TMeasurement, class TFrequencyContainer >
void
Histogram<TMeasurement, TFrequencyContainer >
::SetBinMax( unsigned int dimension, InstanceIdentifier nbin,
                 MeasurementType max)
{
  m_Max[dimension][nbin] = max;
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::BinMinVectorType &
Histogram<TMeasurement, TFrequencyContainer >
::GetDimensionMins(unsigned int dimension) const
{
  return m_Min[dimension];
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::BinMaxVectorType &
Histogram<TMeasurement, TFrequencyContainer >
::GetDimensionMaxs(unsigned int dimension) const
{
  return m_Max[dimension];
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::BinMinContainerType &
Histogram<TMeasurement, TFrequencyContainer >
::GetMins() const
{
  return m_Min;
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram<TMeasurement, TFrequencyContainer >::BinMaxContainerType &
Histogram<TMeasurement, TFrequencyContainer >
::GetMaxs() const
{
  return m_Max;
}

template< class TMeasurement, class TFrequencyContainer >
typename 
Histogram<TMeasurement, TFrequencyContainer >::AbsoluteFrequencyType
Histogram<TMeasurement, TFrequencyContainer >
::GetFrequency( InstanceIdentifier id ) const
{
  return m_FrequencyContainer->GetFrequency(id);
}

template< class TMeasurement, class TFrequencyContainer >
bool 
Histogram<TMeasurement, TFrequencyContainer >
::SetFrequency( InstanceIdentifier id, AbsoluteFrequencyType value)
{
  return m_FrequencyContainer->SetFrequency(id, value);
}

template< class TMeasurement, class TFrequencyContainer >
bool 
Histogram<TMeasurement, TFrequencyContainer >
::IncreaseFrequency(InstanceIdentifier id, AbsoluteFrequencyType value)
{
  return m_FrequencyContainer->IncreaseFrequency(id, value);
}


template< class TMeasurement, class TFrequencyContainer >
void
Histogram<TMeasurement, TFrequencyContainer >
::Initialize(const SizeType &size)
{

  if( this->GetMeasurementVectorSize() == 0 )
    {
    itkExceptionMacro("MeasurementVectorSize is Zero. It should be set to a non-zero value before calling Initialize");
    }

  this->m_Size = size;
  
  // creates offset table which will be used for generation of
  // instance identifiers.
  InstanceIdentifier num = 1;
  
  this->m_OffsetTable.resize( this->GetMeasurementVectorSize() + 1 );

  this->m_OffsetTable[0] = num;
  for (unsigned int i = 0; i < this->GetMeasurementVectorSize(); i++)
    {
    num *= m_Size[i];
    this->m_OffsetTable[i + 1] = num;
    }

  this->m_TempIndex.SetSize( this->GetMeasurementVectorSize() );

  m_NumberOfInstances = num;

  // adjust the sizes of min max value containers 
  unsigned int dim;
  m_Min.resize(this->GetMeasurementVectorSize());
  for ( dim = 0; dim < this->GetMeasurementVectorSize(); dim++)
    {
    m_Min[dim].resize(m_Size[dim]);
    } 

  m_Max.resize(this->GetMeasurementVectorSize());
  for ( dim = 0; dim < this->GetMeasurementVectorSize(); dim++)
    {
    m_Max[dim].resize(m_Size[dim]);
    } 

  // initialize auxiliary variables
  this->m_TempIndex.SetSize( this->GetMeasurementVectorSize() );
  this->m_TempMeasurementVector.SetSize( this->GetMeasurementVectorSize() );

  // initialize the frequency container
  m_FrequencyContainer->Initialize(this->m_OffsetTable[this->GetMeasurementVectorSize()]);
  this->SetToZero();
}

template< class TMeasurement, class TFrequencyContainer >
void 
Histogram<TMeasurement, TFrequencyContainer >
::SetToZero()
{
  m_FrequencyContainer->SetToZero();
}

template< class TMeasurement, class TFrequencyContainer >
void 
Histogram<TMeasurement, TFrequencyContainer >
::Initialize(const SizeType &size, MeasurementVectorType& lowerBound,
             MeasurementVectorType& upperBound) {
  this->Initialize(size);

  float interval;
  for ( unsigned int i = 0; i < this->GetMeasurementVectorSize(); i++)
    {
    if( size[i] > 0 )
      {
      interval = (float) (upperBound[i] - lowerBound[i]) 
                       / static_cast< MeasurementType >(size[i]);

      // Set the min vector and max vector
      for (unsigned int j = 0; j < static_cast<unsigned int>(size[i] - 1); j++)
        {
        this->SetBinMin(i, j, (MeasurementType)(lowerBound[i] +  
                                                ((float)j * interval)));
        this->SetBinMax(i, j, (MeasurementType)(lowerBound[i] +  
                                                (((float)j + 1) * interval)));
        }
      this->SetBinMin(i, size[i] - 1, 
                      (MeasurementType)(lowerBound[i] + 
                                        (((float) size[i] - 1) * interval)));
      this->SetBinMax(i, size[i] - 1, 
                      (MeasurementType)(upperBound[i]));
      }
    }

  //  this->PrintMaxs(std::cout, 3);
  //  this->PrintMins(std::cout, 3);
}

template< class TMeasurement, class TFrequencyContainer >
void 
Histogram<TMeasurement, TFrequencyContainer >
::Initialize(const BinMinContainerType& mins, const BinMaxContainerType& maxs)
{
  if( mins.size() != maxs.size() )
  {
	itkExceptionMacro("MeasurementVectorSize does not match between minimums and maximums vectors");
  }

  if( mins.size() != this->GetMeasurementVectorSize() )
  {
	itkExceptionMacro("MeasurementVectorSize does not match with the one of minimums or maximums vector");
  }

  SizeType size;
  size.SetSize( this->GetMeasurementVectorSize() );
  
  for ( unsigned int i = 0; i< this->GetMeasurementVectorSize(); i++ )
  {
	size[i] = mins[i].size();
  }

  this->Initialize( size );

  m_Min = mins;

  m_Max = maxs;

}

template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram<TMeasurement, TFrequencyContainer >::IndexType &
Histogram<TMeasurement, TFrequencyContainer >
::GetIndex(const MeasurementVectorType& measurement) const
{
  // Have this deprecated method call the un-deprecated one.. 
  this->GetIndex( measurement, m_TempIndex );
  return m_TempIndex;
}


/** */
template< class TMeasurement, class TFrequencyContainer >
bool Histogram<TMeasurement, TFrequencyContainer >
::GetIndex(const MeasurementVectorType & measurement,IndexType & index ) const
{
  // now using something similar to binary search to find
  // index.
  unsigned int dim;

  if( index.Size() != this->GetMeasurementVectorSize() )
    {
    index.SetSize( this->GetMeasurementVectorSize() );
    }
  
  int begin;
  int mid;
  int end;

  MeasurementType median;
  MeasurementType tempMeasurement;

  for (dim = 0; dim < this->GetMeasurementVectorSize(); dim++)
    {
    tempMeasurement = measurement[dim];
    begin = 0;
    if (tempMeasurement < m_Min[dim][begin])
      {
      // one of measurement is below the minimum
      // its ok if we extend the bins to infinity.. not ok if we don't
      if(!m_ClipBinsAtEnds)
        {
        index[dim] = (long) 0;
        continue;
        }
      else
        { // set an illegal value and return 0
        index[dim] = (long) m_Size[dim]; 
        return false;
        }
      }

    end = m_Min[dim].size() - 1;
    if (tempMeasurement >= m_Max[dim][end])
      {
      // one of measurement is above the maximum
      // its ok if we extend the bins to infinity.. not ok if we don't
      if(!m_ClipBinsAtEnds)
        {
        index[dim] = (long) m_Size[dim]-1;
        continue;
        }
      else
        { // set an illegal value and return 0
        index[dim] = (long) m_Size[dim]; 
        return false;
        }
      }

    // Binary search for the bin where this measurement could be
    mid = (end + 1) / 2;
    median = m_Min[dim][mid];

    while(true)
      {
      if (tempMeasurement < median )
        {
        end = mid - 1;
        } 
      else if (tempMeasurement > median)
        {
        // test whether it is inside the current bin by comparing to the max of this bin.
        if( tempMeasurement <  m_Max[dim][mid] && 
            tempMeasurement >= m_Min[dim][mid] )
          {
          index[dim] = mid;
          break;
          }
        // otherwise, continue binary search
        begin = mid + 1;
        }
      else
        {
        index[dim] = mid;
        break;
        }
      mid = begin + (end - begin) / 2;
      median = m_Min[dim][mid];
      } // end of while
    } // end of for()
  return true;
}


template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram<TMeasurement, TFrequencyContainer >::IndexType&
Histogram<TMeasurement, TFrequencyContainer >
::GetIndex( InstanceIdentifier id ) const
{
  InstanceIdentifier id2 = id;

  for (int i = this->GetMeasurementVectorSize() - 1; i > 0; i--)
    {
    m_TempIndex[i] = static_cast<IndexValueType>(id2 / this->m_OffsetTable[i]);
    id2 -= (m_TempIndex[i] * this->m_OffsetTable[i]);
    }
  m_TempIndex[0] = static_cast<IndexValueType>(id2);
  
  return m_TempIndex;
}


template< class TMeasurement, class TFrequencyContainer >
inline bool
Histogram<TMeasurement, TFrequencyContainer >
::IsIndexOutOfBounds(const IndexType &index) const
{
  for (unsigned int dim = 0; dim < this->GetMeasurementVectorSize(); dim++)
    {
    if (index[dim] < 0 || index[dim] >= static_cast<IndexValueType>(m_Size[dim]))
      {
      return true;
      }
    }
  return false;
}

template< class TMeasurement, class TFrequencyContainer >
inline typename Histogram<TMeasurement, TFrequencyContainer >::InstanceIdentifier
Histogram<TMeasurement, TFrequencyContainer >
::GetInstanceIdentifier(const IndexType &index) const
{
  InstanceIdentifier instanceId = 0;
  for (int i= this->GetMeasurementVectorSize() - 1; i > 0; i-- )
    {
    instanceId += index[i] * this->m_OffsetTable[i];
    }
  
  instanceId += index[0];
  
  return instanceId;
}


template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram<TMeasurement, TFrequencyContainer >::MeasurementType&
Histogram<TMeasurement, TFrequencyContainer >
::GetBinMinFromValue(const unsigned int dimension, const float value ) const
{
  // If the value is lower than any of min value in the Histogram,
  // it returns the lowest min value
  if ( value <= this->m_Min[dimension][0] )
    {
    return this->m_Min[dimension][0];
    }

  // If the value is higher than any of min value in the Histogram,
  // it returns the highest min value
  if ( value >= m_Min[dimension][m_Size[dimension]-1] )
    {
    return m_Min[dimension][this->m_Size[dimension]-1];
    }

  unsigned int binMinFromValue = 0;

  for ( unsigned int i = 0; i < this->m_Size[dimension]; i++ )
    {
    if (  (value >= this->m_Min[dimension][i])
          && (value <  this->m_Max[dimension][i])  )
      {
      binMinFromValue = i;
      }
    }

  return this->m_Min[dimension][binMinFromValue];
}

template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram< TMeasurement, TFrequencyContainer >::MeasurementType&
Histogram< TMeasurement, TFrequencyContainer >
::GetBinMaxFromValue(const unsigned int dimension, const float value ) const
{
  // If the value is lower than any of max value in the Histogram,
  // it returns the lowest max value
  if ( value <= this->m_Max[dimension][0] )
    {
    return this->m_Max[dimension][0];
    }

  // If the value is higher than any of max value in the Histogram,
  // it returns the highest max value
  if ( value >= m_Max[dimension][m_Size[dimension]-1] )
    {
    return m_Max[dimension][this->m_Size[dimension]-1];
    }

  unsigned int binMaxFromValue = 0;

  for ( unsigned int i = 0; i < this->m_Size[dimension]; i++ )
    {
    if (  (value >= this->m_Min[dimension][i])
          && (value <  this->m_Max[dimension][i])  )
      {
      binMaxFromValue = i;
      }
    }

  return this->m_Max[dimension][binMaxFromValue];
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram< TMeasurement, TFrequencyContainer >::MeasurementVectorType&
Histogram< TMeasurement, TFrequencyContainer >
::GetHistogramMinFromIndex(const IndexType &index) const
{
  for( unsigned int i=0; i < this->GetMeasurementVectorSize(); i++ )
    {
    m_TempMeasurementVector[i] = this->GetBinMin(i, index[i]);
    }
  return m_TempMeasurementVector;
}

template< class TMeasurement, class TFrequencyContainer >
const typename 
Histogram< TMeasurement, TFrequencyContainer >::MeasurementVectorType&
Histogram< TMeasurement, TFrequencyContainer >
::GetHistogramMaxFromIndex(const IndexType &index) const
{
  for( unsigned int i=0; i < this->GetMeasurementVectorSize(); i++ )
    {
    m_TempMeasurementVector[i] = this->GetBinMax(i, index[i]);
    }
  return m_TempMeasurementVector;
}

template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram< TMeasurement, TFrequencyContainer >::MeasurementVectorType &
Histogram< TMeasurement, TFrequencyContainer >
::GetMeasurementVector(const IndexType &index) const
{
  for ( unsigned int i = 0; i < this->GetMeasurementVectorSize(); i++)
    {
    MeasurementType value = (m_Min[i][index[i]] + m_Max[i][index[i]]);
    m_TempMeasurementVector[i] =  static_cast< MeasurementType >( value / 2.0 );
    }
  return m_TempMeasurementVector;
}

template< class TMeasurement, class TFrequencyContainer >
inline const typename Histogram< TMeasurement, TFrequencyContainer >::MeasurementVectorType &
Histogram< TMeasurement, TFrequencyContainer >
::GetMeasurementVector( InstanceIdentifier id) const
{
  return this->GetMeasurementVector( this->GetIndex(id) );
}

template< class TMeasurement, class TFrequencyContainer >
inline void
Histogram< TMeasurement, TFrequencyContainer >
::SetFrequency(const AbsoluteFrequencyType value) 
{
  typename Self::Iterator iter = this->Begin();
  typename Self::Iterator end = this->End();
  
  while ( iter != end )
    {
    iter.SetFrequency(value);
    ++iter;
    }
}

template< class TMeasurement, class TFrequencyContainer >
inline bool
Histogram< TMeasurement, TFrequencyContainer >
::SetFrequency(const IndexType &index, const AbsoluteFrequencyType value) 
{
  return this->SetFrequency( this->GetInstanceIdentifier(index), value);
}
  
template< class TMeasurement, class TFrequencyContainer >
inline bool
Histogram< TMeasurement, TFrequencyContainer >
::SetFrequency(const MeasurementVectorType &measurement, const AbsoluteFrequencyType value) 
{
  return this->SetFrequency( this->GetInstanceIdentifier(GetIndex(measurement)), value);
}

template< class TMeasurement, class TFrequencyContainer >
inline bool
Histogram< TMeasurement, TFrequencyContainer >
::IncreaseFrequency(const IndexType &index, const AbsoluteFrequencyType value)
{
  const bool result = 
      this->IncreaseFrequency( this->GetInstanceIdentifier(index), value);
  return result;
}
  
template< class TMeasurement, class TFrequencyContainer >
inline bool
Histogram< TMeasurement, TFrequencyContainer >
::IncreaseFrequency(const MeasurementVectorType &measurement, const AbsoluteFrequencyType value) 
{
  IndexType index;
  this->GetIndex( measurement, index );
  return this->IncreaseFrequency( this->GetInstanceIdentifier( index ), value );
}


template< class TMeasurement, class TFrequencyContainer >
inline typename Histogram< TMeasurement, TFrequencyContainer >::AbsoluteFrequencyType
Histogram< TMeasurement, TFrequencyContainer >
::GetFrequency(const IndexType &index) const
{
  return ( this->GetFrequency( this->GetInstanceIdentifier(index)) );
}

template< class TMeasurement, class TFrequencyContainer >
typename 
Histogram< TMeasurement, TFrequencyContainer >::MeasurementType 
Histogram< TMeasurement, TFrequencyContainer >
::GetMeasurement( InstanceIdentifier n, unsigned int dimension) const
{
  return static_cast< MeasurementType >((m_Min[dimension][n] + 
                                         m_Max[dimension][n]) / 2); 
}

template< class TMeasurement, class TFrequencyContainer >
typename 
Histogram< TMeasurement, TFrequencyContainer >::AbsoluteFrequencyType
Histogram< TMeasurement, TFrequencyContainer >
::GetFrequency( InstanceIdentifier n, unsigned int dimension) const
{
  InstanceIdentifier nextOffset = this->m_OffsetTable[dimension + 1];
  InstanceIdentifier current = this->m_OffsetTable[dimension] * n;
  InstanceIdentifier includeLength = this->m_OffsetTable[dimension];
  InstanceIdentifier include;
  InstanceIdentifier includeEnd;
  InstanceIdentifier last = this->m_OffsetTable[this->GetMeasurementVectorSize()];

  AbsoluteFrequencyType frequency = 0;
  while (current < last)
    {
    include = current;
    includeEnd = include + includeLength;
    while(include < includeEnd)
      {
      frequency += GetFrequency(include);
      include++;
      }
    current += nextOffset;
    }
  return frequency;
}

template< class TMeasurement, class TFrequencyContainer >
inline typename Histogram< TMeasurement, TFrequencyContainer >::TotalAbsoluteFrequencyType
Histogram< TMeasurement, TFrequencyContainer >
::GetTotalFrequency() const
{
  return m_FrequencyContainer->GetTotalFrequency();
}

template< class TMeasurement, class TFrequencyContainer >
double
Histogram< TMeasurement, TFrequencyContainer >
::Quantile(unsigned int dimension, double p) const
{
  InstanceIdentifier n;
  const unsigned int size = this->GetSize(dimension);
  double p_n_prev;
  double p_n;
  double f_n;
  double cumulated = 0;
  double totalFrequency = double( this->GetTotalFrequency() );
  double binProportion;
  double min, max, interval;

  if ( p < 0.5 )
    {
    n = 0;
    p_n = NumericTraits< double >::Zero;
    do 
      {
      f_n = this->GetFrequency(n, dimension);
      cumulated += f_n;
      p_n_prev = p_n;
      p_n = cumulated / totalFrequency;
      n++;
      } 
    while( n < size && p_n < p);

    binProportion = f_n / totalFrequency;

    min = double( this->GetBinMin(dimension, n - 1) );
    max = double( this->GetBinMax(dimension, n - 1) );
    interval = max - min;
    return min + ((p - p_n_prev) / binProportion) * interval;
    }
  else
    {
    n = size - 1;
    InstanceIdentifier m = NumericTraits< InstanceIdentifier >::Zero;
    p_n      = NumericTraits< double >::One;
    do 
      {
      f_n = this->GetFrequency(n, dimension);
      cumulated += f_n;
      p_n_prev = p_n;
      p_n = NumericTraits< double >::One - cumulated / totalFrequency;
      n--;
      m++;
      } 
    while( m < size && p_n > p);

    binProportion = f_n / totalFrequency;
    min = double( this->GetBinMin(dimension, n + 1) );
    max = double( this->GetBinMax(dimension, n + 1) );
    interval = max - min;
    return max - ((p_n_prev - p) / binProportion) * interval;
    }
}

template< class TMeasurement, class TFrequencyContainer >
void 
Histogram< TMeasurement, TFrequencyContainer >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "OffsetTable: " <<  std::endl;
  for(unsigned int i=0; i < this->m_OffsetTable.size(); i++)
    {
    os << this->m_OffsetTable[i] << "  ";
    }
  os << std::endl;
  os << indent << "ClipBinsAtEnds: " <<  
    itk::NumericTraits<bool>::PrintType( this->GetClipBinsAtEnds() ) << std::endl;
  os << indent << "FrequencyContainerPointer: " << m_FrequencyContainer
     << std::endl;
}

template< class TMeasurement, class TFrequencyContainer >
void 
Histogram< TMeasurement, TFrequencyContainer >
::PrintMaxs(std::ostream& os, Indent indent) const
{
	os << indent << "Maxs:" << std::endl;

	for(unsigned int i=0; i<  this->GetMeasurementVectorSize(); i++) {
		  os << indent << "Dimension " << i << ":" << std::endl;
		  BinMaxVectorType maxs = this->GetDimensionMaxs(i);

		  os << indent << indent << "Maxs: [ ";
		  for ( unsigned int j=0; j< maxs.size(); j++)  os << maxs[j] << " ";
		  os << "]" << std::endl;
	  }
}

template< class TMeasurement, class TFrequencyContainer >
void
Histogram< TMeasurement, TFrequencyContainer >
::PrintMins(std::ostream& os, Indent indent) const
{
	os << indent << "Mins:" << std::endl;

	for(unsigned int i=0; i<  this->GetMeasurementVectorSize(); i++) {
		  os << indent << "Dimension " << i << ":" << std::endl;
		  BinMinVectorType mins = this->GetDimensionMins(i);

		  os << indent << indent << "Mins: [ ";
		  for ( unsigned int j=0; j< mins.size(); j++)  os << mins[j] << " ";
		  os << "]" << std::endl;
	  }
}
template< class TMeasurement, class TFrequencyContainer >
void
Histogram< TMeasurement, TFrequencyContainer >
::Graft( const DataObject *thatObject )
{
  this->Superclass::Graft(thatObject);

  const Self *thatConst = dynamic_cast< const Self * >(thatObject);
  if (thatConst)
    {
    Self *that = const_cast< Self * >(thatConst); 
    this->m_Size                  = that->m_Size;
    this->m_OffsetTable           = that->m_OffsetTable;
    this->m_FrequencyContainer    = that->m_FrequencyContainer;
    this->m_NumberOfInstances     = that->m_NumberOfInstances;
    this->m_Min                   = that->m_Min;
    this->m_Max                   = that->m_Max;
    this->m_TempMeasurementVector = that->m_TempMeasurementVector;
    this->m_TempIndex             = that->m_TempIndex;
    this->m_ClipBinsAtEnds        = that->m_ClipBinsAtEnds;
    }
}

} // end of namespace Statistics 
} // end of namespace itk 

#endif
