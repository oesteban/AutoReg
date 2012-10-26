/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSampleToAnisotropicHistogramFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009-05-02 05:43:58 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSampleToAnisotropicHistogramFilter_txx
#define __itkSampleToAnisotropicHistogramFilter_txx

#include "itkSampleToAnisotropicHistogramFilter.h"
#include "itkStatisticsAlgorithm.h"

namespace itk { 
namespace Statistics {

template < class TSample, class THistogram >
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::SampleToAnisotropicHistogramFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs( 1 );
  this->ProcessObject::SetNumberOfRequiredOutputs(1);

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0) );

  const unsigned int minimumNumberOfComponents = 1;

  // Set some default inputs 
  HistogramSizeType histogramSize( minimumNumberOfComponents );
  histogramSize.Fill(0);
  this->SetHistogramSize( histogramSize );

  this->SetMarginalScale( 100 );

  this->SetAutoMinimumMaximum( true );

  this->SetPDFDensityFactor( 100 );
  this->SetMaxEstimatorIterations( 10 );
}

template < class TSample, class THistogram >
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::~SampleToAnisotropicHistogramFilter()
{
}

template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::SetInput( const SampleType * sample )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, 
                                   const_cast< SampleType * >( sample ) );
}

template < class TSample, class THistogram >
const typename 
SampleToAnisotropicHistogramFilter< TSample, THistogram >::SampleType *
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GetInput() const
{
  const SampleType * input = 
    static_cast< const SampleType * >( this->ProcessObject::GetInput( 0 ) );
  return input;
}


template < class TSample, class THistogram >
const typename 
SampleToAnisotropicHistogramFilter< TSample, THistogram >::HistogramType *
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GetOutput() const
{
  const HistogramType * output = 
    static_cast<const HistogramType*>(this->ProcessObject::GetOutput(0));

  return output;
}

template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GraftOutput(DataObject *graft)
{
  DataObject * output = 
   const_cast< HistogramType * >( this->GetOutput() );

  // Call Histogram to copy meta-information, and the container
  output->Graft( graft );
}


template < class TSample, class THistogram >
typename SampleToAnisotropicHistogramFilter< TSample, THistogram >::DataObjectPointer
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::MakeOutput(unsigned int)
{
  return static_cast<DataObject*>(HistogramType::New().GetPointer());
}


template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  // m_AutoMinimumMaximum
  os << indent << "AutoMinimumMaximum: "
     << this->GetAutoMinimumMaximumInput() << std::endl;
  // m_MarginalScale
  os << indent << "MarginalScale: "
     << this->GetMarginalScaleInput() << std::endl;
  // outputHistogramBinMinimum
  os << indent << "HistogramBinMinimum: "
     << this->GetHistogramBinMinimumInput() << std::endl;
  // outputHistogramBinMaximum
  os << indent << "HistogramBinMaximum: "
     << this->GetHistogramBinMaximumInput() << std::endl;
  // outputHistogramSize
  os << indent << "HistogramSize: "
     << this->GetHistogramSizeInput() << std::endl;
}


template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GeneratePDF(HistogramMeasurementVectorType lower, HistogramMeasurementVectorType upper)
{
	pdfGenerator = PDFHistogramFilter::New();
	pdfGenerator->SetInput( this->GetInput() );

	HistogramSizeType size = this->GetHistogramSizeInput()->Get();
	unsigned int factor = this->GetPDFDensityFactorInput()->Get();
	for ( unsigned int i = 0; i < size.size() ; i++)
		size[i] = size[i] * factor;


	pdfGenerator->SetHistogramSize( size );
	pdfGenerator->SetAutoMinimumMaximumInput( false );

	pdfGenerator->SetMarginalScaleInput( this->GetMarginalScaleInput() );

	pdfGenerator->SetHistogramBinMinimum( lower );
	pdfGenerator->SetHistogramBinMaximum( upper );
	pdfGenerator->Update();
	m_pdf = pdfGenerator->GetOutput();
}

template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::ComputeBinWidths( unsigned int dim )
{
	HistogramType * outputHistogram = 
		static_cast<HistogramType*>(this->ProcessObject::GetOutput(0));

	typename HistogramType::ConstIterator pdf_it = m_pdf->Begin();

	typename HistogramType::Iterator hist_it = outputHistogram->Begin();
	typename HistogramType::Iterator hist_end = outputHistogram->End();

	typename HistogramType::BinMaxVectorType maxs = outputHistogram->GetDimensionMaxs( dim );
	typename HistogramType::BinMinVectorType mins = outputHistogram->GetDimensionMins( dim );

	// Calculate centroids
	std::vector< HistogramMeasurementType > centroids;
	while ( hist_it != hist_end )
	{
		typename HistogramType::InstanceIdentifier id = hist_it.GetInstanceIdentifier();

		HistogramMeasurementType binWeightedSum=0;
		HistogramMeasurementType binSum=0;

		while( pdf_it != m_pdf->End() &&
			m_pdf->GetBinMax( dim , pdf_it.GetInstanceIdentifier() ) <= maxs[id] )
		{
			binSum+= pdf_it.GetFrequency();
			HistogramMeasurementVectorType value = pdf_it.GetMeasurementVector();

			binWeightedSum += pdf_it.GetFrequency() * value[dim];
			++pdf_it;
		}

		HistogramMeasurementType centroid = (binSum>0)?(binWeightedSum/binSum):((maxs[id]+mins[id])*0.5);
		centroids.push_back( centroid );

		++hist_it;
	}

	// Bin Maxs/Mins update
	hist_it = outputHistogram->Begin();
	unsigned int end_id = hist_end.GetInstanceIdentifier();
	while ( hist_it != hist_end )
	{
		unsigned int id = hist_it.GetInstanceIdentifier();

		if (id > 0) {
			outputHistogram->SetBinMin( dim, id, outputHistogram->GetBinMax( dim, id-1 ) );
		}

		HistogramMeasurementType upperBound;

		if ( id < end_id-1 ){
			upperBound = (centroids[id]+centroids[id+1])*0.5;
		}
		else 
		{
			upperBound = outputHistogram->GetBinMax( dim, id );
		}

		if ( upperBound <= outputHistogram->GetBinMin( dim, id ) )
		{
			itkExceptionMacro( << "Trying to set maximum of bin id=" << id << " at dim=" <<dim << " lower than its minimum" );
		}

		outputHistogram->SetBinMax( dim, id, upperBound );

		++hist_it;
	}
}

template < class TSample, class THistogram >
double
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GetNoisePower( unsigned int dim )
{
	HistogramType * outputHistogram = 
		static_cast<HistogramType*>(this->ProcessObject::GetOutput(0));

	typename HistogramType::ConstIterator pdf_it = m_pdf->Begin();
	typename HistogramType::ConstIterator pdf_end = m_pdf->End();

	typename HistogramType::Iterator hist_it = outputHistogram->Begin();
	typename HistogramType::Iterator hist_end = outputHistogram->End();

	typename HistogramType::BinMaxVectorType maxs = outputHistogram->GetDimensionMaxs( dim );

	double noisePowerSum=0;

	while ( hist_it != hist_end )
	{
		typename HistogramType::InstanceIdentifier id = hist_it.GetInstanceIdentifier();
		HistogramMeasurementType c = outputHistogram->GetMeasurement(id, dim);

		while( pdf_it != pdf_end &&
			m_pdf->GetBinMax( dim , pdf_it.GetInstanceIdentifier() ) <= maxs[id] )
		{
			HistogramMeasurementType x = m_pdf->GetMeasurement(pdf_it.GetInstanceIdentifier(), dim);
			noisePowerSum+= pow( (c - x), 2) * pdf_it.GetFrequency();
			++pdf_it;
		}
		++hist_it;
	}
	noisePowerSum /= m_pdf->GetTotalFrequency();

	return noisePowerSum;

}

template < class TSample, class THistogram >
double
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GetSignalPower( unsigned int dim )
{
	double signalPower = 0;

	typename HistogramType::ConstIterator pdf_it = m_pdf->Begin();
	typename HistogramType::ConstIterator pdf_end = m_pdf->End();

	while( pdf_it != pdf_end )
	{
		HistogramMeasurementType x = m_pdf->GetMeasurement(pdf_it.GetInstanceIdentifier(), dim);
		signalPower+= pow( x, 2) * pdf_it.GetFrequency();
		++pdf_it;
	}

	signalPower /= m_pdf->GetTotalFrequency();

	return signalPower;
}


template < class TSample, class THistogram >
double
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GetErrorSum( unsigned int dim )
{

	HistogramType * outputHistogram = 
		static_cast<HistogramType*>(this->ProcessObject::GetOutput(0));
	
	unsigned int hTotalFreq = outputHistogram->GetTotalFrequency();

	if ( hTotalFreq == 0 ){
		this->Sample();
		outputHistogram->GetTotalFrequency();
	}
	

	typename HistogramType::ConstIterator pdf_it = m_pdf->Begin();
	typename HistogramType::ConstIterator pdf_end = m_pdf->End();

	typename HistogramType::Iterator hist_it = outputHistogram->Begin();
	typename HistogramType::Iterator hist_end = outputHistogram->End();

	typename HistogramType::BinMaxVectorType pdf_maxs = m_pdf->GetDimensionMaxs( dim );
	typename HistogramType::BinMinVectorType pdf_mins = m_pdf->GetDimensionMins( dim );

	typename HistogramType::BinMaxVectorType maxs = outputHistogram->GetDimensionMaxs( dim );
	typename HistogramType::BinMinVectorType mins = outputHistogram->GetDimensionMins( dim );

	double errorSum=0;

	
	unsigned int pdfTotalFreq = m_pdf->GetTotalFrequency();

	double totalPx=0;
	double totalPh=0;

	while ( hist_it != hist_end )
	{
		typename HistogramType::InstanceIdentifier id = hist_it.GetInstanceIdentifier();
		HistogramMeasurementType c = outputHistogram->GetMeasurement(id, dim);


		double pdfBinSum = 0;
		HistogramMeasurementType Wpdf = 0;
		unsigned int pdfBins = 0;
		double binError = 0;
		double binTotalPx=0;
		double binTotalPh= 0;

		unsigned int binTotalFx= 0;
		unsigned int binTotalFh= 0;


		unsigned int Fh = hist_it.GetFrequency();
		//unsigned int Fhx= Fh/100;

		double Ph = (1.0 * Fh) / hTotalFreq;
		//double Phx = Ph/100;
		HistogramMeasurementType teoric_wh = maxs[id]-mins[id];
		HistogramMeasurementType wh = 0;


		while( pdf_it != pdf_end &&
			pdf_maxs[pdf_it.GetInstanceIdentifier()] <= maxs[id] )
		{
			unsigned int pdf_id = pdf_it.GetInstanceIdentifier();
			HistogramMeasurementType x = m_pdf->GetMeasurement(pdf_id, dim);
			HistogramMeasurementType w = pdf_maxs[pdf_id] - pdf_mins[pdf_id];
			wh+=w;
			unsigned int Fx = pdf_it.GetFrequency();
			double Px = (1.0*Fx)/pdfTotalFreq;
			Wpdf += w;
			pdfBinSum+= w * Px;
			pdfBins++;
			//binError+= (Px - Phx) * w;
			binTotalPx+= Px;
			//binTotalPh+=Phx;
			++pdf_it;



			//					binTotalFh+=Fhx;
			//					binTotalFx+= Fx;
		}

		totalPx+= binTotalPx;
		totalPh+=Ph;

		double sampleError = binTotalPx - Ph;
		binError = sampleError * (wh/pdfBins);
		errorSum += binError;

		//std::cout << "Error(" << id << ")=" << binError << ". binPx="<< binTotalPx << ", Ph=" << Ph << ", Phx="<< binTotalPh << ", sampleError=" << sampleError <<", wh=" <<wh << "/" << teoric_wh << std::endl;
		////std::cout << "Error(" << id << ")=" << binError << ". binFx="<< binTotalFx << ", Fh=" << Fh << ", Fhx="<< binTotalFh << ", bins=" <<pdfBins << std::endl;
		//if ( abs(wh - Wpdf)> Wpdf/pdfBins )
		//	std::cout << "Width mismatch Wh=" << wh << ", Wpdf=" << Wpdf << "." << std::endl;
		///*std::cout << "Frequencies differ at bin " << id << "(" << pdfSum - Ph << ")" << std::endl;*/
		++hist_it;
	}

	//std::cout << "Total error=" << errorSum << std::endl;

	return errorSum;
}

template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::Sample()
{
  const SampleType * inputSample = this->GetInput();
  const typename SampleType::InstanceIdentifier measurementVectorSize =
    inputSample->GetMeasurementVectorSize();
  HistogramType * outputHistogram = 
    static_cast<HistogramType*>(this->ProcessObject::GetOutput(0));

  outputHistogram->SetToZero();

  typename SampleType::ConstIterator iter = inputSample->Begin();
  typename SampleType::ConstIterator last = inputSample->End();
  
  typename SampleType::MeasurementVectorType lvector;
  
  typename HistogramType::IndexType index( measurementVectorSize);  
  typename HistogramType::MeasurementVectorType hvector( measurementVectorSize );
  
  unsigned int i;
  
  while (iter != last)
    {
    lvector = iter.GetMeasurementVector();
    for ( i = 0; i < inputSample->GetMeasurementVectorSize(); i++)
      {
      hvector[i] = static_cast< HistogramMeasurementType >( lvector[i] );
      }

    outputHistogram->GetIndex( hvector, index );
    if( !outputHistogram->IsIndexOutOfBounds( index ) )
      {
      // if the measurement vector is out of bound then
      // the GetIndex method has returned an index set to the max size of
      // the invalid dimension - even if the hvector is less than the minimum
      // bin value.
      // If the index isn't valid, we don't increase the frequency.
      // See the comments in Histogram->GetIndex() for more info.
      outputHistogram->IncreaseFrequency( index, 1 );
      }
    ++iter;
    }
}

template < class TSample, class THistogram >
void
SampleToAnisotropicHistogramFilter< TSample, THistogram >
::GenerateData()
{
  const SampleType * inputSample = this->GetInput();

  const InputHistogramMeasurementVectorObjectType * binMinimumObject =
    this->GetHistogramBinMinimumInput();

  const InputHistogramMeasurementVectorObjectType * binMaximumObject =
    this->GetHistogramBinMaximumInput();

  const InputHistogramMeasurementObjectType * marginalScaleObject =
    this->GetMarginalScaleInput();
 
  const InputBooleanObjectType * autoMinimumMaximum =
    this->GetAutoMinimumMaximumInput();

  const InputHistogramSizeObjectType * histogramSizeObject =
    this->GetHistogramSizeInput();
    
  if( histogramSizeObject == NULL )
    {
    itkExceptionMacro("Histogram Size input is missing");
    }
    
  if( marginalScaleObject == NULL )
    {
    itkExceptionMacro("Marginal scale input is missing");
    }

  HistogramSizeType histogramSize = histogramSizeObject->Get();

  HistogramMeasurementType marginalScale = marginalScaleObject->Get();

  HistogramType * outputHistogram = 
    static_cast<HistogramType*>(this->ProcessObject::GetOutput(0));

  const typename SampleType::InstanceIdentifier measurementVectorSize =
    inputSample->GetMeasurementVectorSize();

  if( measurementVectorSize == 0 )
    {
    itkExceptionMacro("Input sample MeasurementVectorSize is zero");
    }

  if( histogramSize.Size() != measurementVectorSize )
    {
    itkExceptionMacro("Histogram number of components: "
                      << histogramSize.Size()
                      << " doesn't match Measurement Vector Size: "
                      << measurementVectorSize);
    }

  outputHistogram->SetMeasurementVectorSize( measurementVectorSize );

  typename SampleType::MeasurementVectorType lower;
  typename SampleType::MeasurementVectorType upper;

  MeasurementVectorTraits::SetLength( lower, measurementVectorSize );
  MeasurementVectorTraits::SetLength( upper, measurementVectorSize );
 
  HistogramMeasurementVectorType h_upper;
  HistogramMeasurementVectorType h_lower;

  MeasurementVectorTraits::SetLength( h_lower, measurementVectorSize );
  MeasurementVectorTraits::SetLength( h_upper, measurementVectorSize );

  const HistogramMeasurementType maximumPossibleValue = 
    itk::NumericTraits< HistogramMeasurementType >::max();

  if( autoMinimumMaximum && autoMinimumMaximum->Get() )
    {
    if( inputSample->Size() )
      {
      Algorithm::FindSampleBound( 
        inputSample,  inputSample->Begin(), inputSample->End(), lower, upper);
      
      for( unsigned int i = 0; i < measurementVectorSize; i++ )
        {
        if( !NumericTraits< HistogramMeasurementType >::is_integer )
          {
          const double margin = 
              ( static_cast< HistogramMeasurementType >( upper[i] - lower[i] ) / 
                static_cast< HistogramMeasurementType >( histogramSize[i] ) ) / 
              static_cast< HistogramMeasurementType >( marginalScale );

          // Now we check if the upper[i] value can be increased by 
          // the margin value without saturating the capacity of the 
          // HistogramMeasurementType
          if( ( maximumPossibleValue - upper[i] ) > margin )
            {
            h_upper[i] = static_cast< HistogramMeasurementType > (upper[i] + margin);
            }
          else
            { 
            // an overflow would occur if we add 'margin' to the upper 
            // therefore we just compromise in setting h_upper = upper.
            h_upper[i] = static_cast< HistogramMeasurementType >( upper[i] );
            // Histogram measurement type would force the clipping the max value.
            // Therefore we must call the following to include the max value:
            outputHistogram->SetClipBinsAtEnds(false);
            // The above function is okay since here we are within the autoMinMax 
            // computation and clearly the user intended to include min and max.
            }
          }
        else
          {
          h_upper[i] = (static_cast< HistogramMeasurementType >( upper[i] ) ) + 
            NumericTraits< HistogramMeasurementType  >::One;

          if( h_upper[i] <= upper[i] )
            { 
            // an overflow has occurred therefore set upper to upper
            h_upper[i] = static_cast< HistogramMeasurementType >( upper[i] );
            // Histogram measurement type would force the clipping the max value.
            // Therefore we must call the following to include the max value:
            outputHistogram->SetClipBinsAtEnds(false);
            // The above function is okay since here we are within the autoMinMax 
            // computation and clearly the user intended to include min and max.
            }
          }
        h_lower[i] = static_cast< HistogramMeasurementType >( lower[i] );
        }
      }
    else
      {
      for( unsigned int i = 0; i < measurementVectorSize; i++ )
        {
        h_lower[i] = static_cast< HistogramMeasurementType >( lower[i] );
        h_upper[i] = static_cast< HistogramMeasurementType >( upper[i] );
        }
      }
    }
  else
    {
    if( binMaximumObject == NULL )
      {
      itkExceptionMacro("Histogram Bin Maximum input is missing");
      }

    if( binMinimumObject == NULL )
      {
      itkExceptionMacro("Histogram Bin Minimum input is missing");
      }

    h_upper = binMaximumObject->Get();
    h_lower = binMinimumObject->Get();
    }

  // PDF Generation
  this->GeneratePDF( h_lower, h_upper );

  // initialize the Histogram object using the sizes and
  // the upper and lower bound from the FindSampleBound function
  outputHistogram->Initialize( histogramSize, h_lower, h_upper );

  for ( unsigned int i = 0; i < measurementVectorSize; i++ )
  {
	  unsigned int optEndCondition = this->GetMaxEstimatorIterationsInput()->Get();
	  // Iterative calculation of new bin sizes
	  while( optEndCondition )
	  {
		this->ComputeBinWidths( i );
		optEndCondition--;
		// TODO: noise power condition
	  }
  }

  this->Sample();
}

} // end of namespace Statistics 
} // end of namespace itk

#endif
