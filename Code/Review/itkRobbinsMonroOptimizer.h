/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRobbinsMonroOptimizer.h,v $
  Language:  C++
  Date:      $Date: 2009-06-24 12:02:51 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRobbinsMonroOptimizer_h
#define __itkRobbinsMonroOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"
#include <string>
namespace itk
{
  
/** \class RobbinsMonroOptimizer
 * \brief Implement a Robbins Monro gradient descent optimizer
 *
 * RobbinsMonroOptimizer implements a simple gradient descent optimizer.
 * At each iteration the current position is updated according to
 *
 * \f[ 
 *        p_{n+1} = p_n 
 *                + \mbox{gamma_k} 
                  \, \frac{\partial f(p_n) }{\partial p_n} 
                  
          gamma_k = a / (k + A )^alpha
 * \f]
 *
 * The gamma_k parameter is a monotone decreasing function of iteration k.
 * Using members Seta(), SetAlpha() and SetA() then gamma_k is set up.
 * By default, alpha=1 (convergence if k->infinitum).
 *
 * The optimizer steps through a user defined number of iterations;
 * no convergence checking is done.
 *
 * Additionally, user can scale each component of the df / dp
 * but setting a scaling vector using method SetScale().
 *
 * \sa RegularStepRobbinsMonroOptimizer
 * 
 * \ingroup Numerics Optimizers
 */  
class ITK_EXPORT RobbinsMonroOptimizer : 
    public SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef RobbinsMonroOptimizer          Self;
  typedef SingleValuedNonLinearOptimizer    Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( RobbinsMonroOptimizer, SingleValuedNonLinearOptimizer );


  /** Codes of stopping conditions */
  typedef enum
    {
    MaximumNumberOfIterations,
    MetricError
    } StopConditionType;

  /** Methods to configure the cost function. */
  itkGetConstReferenceMacro( Maximize, bool );
  itkSetMacro( Maximize, bool );
  itkBooleanMacro( Maximize );
  bool GetMinimize( ) const
    { return !m_Maximize; }
  void SetMinimize(bool v)
    { this->SetMaximize(!v); }
  void MinimizeOn()
    { this->MaximizeOff(); }
  void MinimizeOff()
    { this->MaximizeOn(); }
  
  /** Advance one step following the gradient direction. */
  virtual void AdvanceOneStep( void );

  /** Start optimization. */
  void    StartOptimization( void );

  /** Resume previously stopped optimization with current parameters
   * \sa StopOptimization. */
  void    ResumeOptimization( void );

  /** Stop optimization.
   * \sa ResumeOptimization */
  void    StopOptimization( void );

  /** Set the a, A and Alpha parameters. */
  itkSetMacro( Alpha, double );
  itkSetMacro ( a, double );
  itkSetMacro ( A, double );
  itkGetConstReferenceMacro( Alpha, double);
  itkGetConstReferenceMacro( a, double );
  itkGetConstReferenceMacro( A, double );

  /** Set the number of iterations. */
  itkSetMacro( NumberOfIterations, unsigned long );

  /** Get the number of iterations. */
  itkGetConstReferenceMacro( NumberOfIterations, unsigned long );

  /** Get the current iteration number. */
  itkGetConstMacro( CurrentIteration, unsigned long );

  /** Get the current value. */
  itkGetConstReferenceMacro( Value, double );

  /** Get Stop condition. */
  itkGetConstReferenceMacro( StopCondition, StopConditionType );
  const std::string GetStopConditionDescription() const;

  /** Get Gradient condition. */
  itkGetConstReferenceMacro( Gradient, DerivativeType );
  
  double GetLearningRate();

protected:
  RobbinsMonroOptimizer();
  virtual ~RobbinsMonroOptimizer() {};
  void PrintSelf(std::ostream& os, Indent indent) const;


  // made protected so subclass can access
  DerivativeType                m_Gradient; 
  bool                          m_Maximize;
  double                        m_Alpha;
  double                        m_a;
  double                        m_A;

private:
  RobbinsMonroOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  bool                          m_Stop;
  double                        m_Value;
  StopConditionType             m_StopCondition;
  unsigned long                 m_NumberOfIterations;
  unsigned long                 m_CurrentIteration;
  OStringStream                 m_StopConditionDescription;
};

} // end namespace itk


#endif
