/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRobbinsMonroOptimizer.cxx,v $
  Language:  C++
  Date:      $Date: 2009-06-24 12:02:51 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkRobbinsMonroOptimizer_txx
#define _itkRobbinsMonroOptimizer_txx

#include "itkRobbinsMonroOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkExceptionObject.h"

namespace itk
{

/**
 * Constructor
 */
RobbinsMonroOptimizer
::RobbinsMonroOptimizer()
{
  itkDebugMacro("Constructor");

  m_Alpha = 1.0;
  m_a =     100.0;
  m_A =     100.0;
  
  m_NumberOfIterations = 100;
  m_CurrentIteration = 0;
  m_Maximize = false;
  m_Value = 0.0;
  m_StopCondition = MaximumNumberOfIterations;
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
}

const std::string
RobbinsMonroOptimizer
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

void
RobbinsMonroOptimizer
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "A: "
     << m_A << std::endl;
  os << indent << "a: "
     << m_a << std::endl;
  os << indent << "Alpha: "
     << m_Alpha << std::endl;
  os << indent << "NumberOfIterations: "
     << m_NumberOfIterations << std::endl;
  os << indent << "Maximize: "
     << m_Maximize << std::endl;
  os << indent << "CurrentIteration: "
     << m_CurrentIteration;
  os << indent << "Value: "
     << m_Value;
  if (m_CostFunction)
    {
    os << indent << "CostFunction: "
       << m_CostFunction;
    }
  os << indent << "StopCondition: "
     << m_StopCondition;
  os << std::endl;
  os << indent << "Gradient: "
     << m_Gradient;
  os << std::endl;
}


/**
 * Start the optimization
 */
void
RobbinsMonroOptimizer
::StartOptimization( void )
{

  itkDebugMacro("StartOptimization");
   
  m_CurrentIteration   = 0;

  this->SetCurrentPosition( this->GetInitialPosition() );
  this->ResumeOptimization();

}

/**
 * Resume the optimization
 */
void
RobbinsMonroOptimizer
::ResumeOptimization( void )
{
  
  itkDebugMacro("ResumeOptimization");

  m_Stop = false;

  m_StopConditionDescription.str("");
  m_StopConditionDescription << this->GetNameOfClass() << ": ";
  InvokeEvent( StartEvent() );
  while( !m_Stop ) 
    {

    try
      {
      m_CostFunction->GetValueAndDerivative( 
        this->GetCurrentPosition(), m_Value, m_Gradient );
      }
    catch( ExceptionObject& err )
      {
      // An exception has occurred. 
      // Terminate immediately.
      m_StopCondition = MetricError;
      m_StopConditionDescription << "Metric error";
      StopOptimization();

      // Pass exception to caller
      throw err;
      }


    if( m_Stop )
      {
      m_StopConditionDescription << "StopOptimization() called";
      break;
      }
  
    AdvanceOneStep();

    m_CurrentIteration++;

    if( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopConditionDescription << "Maximum number of iterations ("
                                 << m_NumberOfIterations
                                 << ") exceeded.";
      m_StopCondition = MaximumNumberOfIterations;
      StopOptimization();
      break;
      }
    }
}

/**
 * Stop optimization
 */
void
RobbinsMonroOptimizer
::StopOptimization( void )
{

  itkDebugMacro("StopOptimization");

  m_Stop = true;
  InvokeEvent( EndEvent() );
}

/**
 * Advance one Step following the gradient direction
 */
void
RobbinsMonroOptimizer
::AdvanceOneStep( void )
{ 

  itkDebugMacro("AdvanceOneStep");

  double direction;
  if( this->m_Maximize ) 
    {
    direction = 1.0;
    }
  else 
    {
    direction = -1.0;
    }

  const unsigned int spaceDimension =  m_CostFunction->GetNumberOfParameters();

  const ParametersType & currentPosition = this->GetCurrentPosition();

  ScalesType scales = this->GetScales();

  // Make sure the scales have been set properly
  if (scales.size() != spaceDimension)
    {
    itkExceptionMacro(<< "The size of Scales is "
                      << scales.size()
                      << ", but the NumberOfParameters for the CostFunction is "
                      << spaceDimension
                      << ".");
    }

  DerivativeType transformedGradient( spaceDimension ); 

  for(unsigned int j = 0; j < spaceDimension; j++)
    {
    transformedGradient[j] = m_Gradient[j] / scales[j];
    }

  ParametersType newPosition( spaceDimension );
  for(unsigned int j = 0; j < spaceDimension; j++)
    {
    newPosition[j] = currentPosition[j] + 
      direction * this->GetLearningRate() * transformedGradient[j];
    }

  this->SetCurrentPosition( newPosition );

  this->InvokeEvent( IterationEvent() );

}

double
RobbinsMonroOptimizer
::GetLearningRate( void )
{
  return (m_a / std::pow( (double) (this->GetCurrentIteration() + m_A), m_Alpha ) );
}

} // end namespace itk

#endif
