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

#ifndef OPTIMIZERHELPER_H
#define OPTIMIZERHELPER_H

#include <boost/any.hpp>
#include <boost/type_traits.hpp>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <iomanip>


#include <itkObject.h>
#include <itkSmartPointer.h>

#include <itkSingleValuedNonLinearOptimizer.h>
#include <itkExhaustiveOptimizer.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkRobbinsMonroOptimizer.h>
#include <itkSPSAOptimizer.h>
#include <itkPowellOptimizer.h>
#include <itkFRPROptimizer.h>
#include <itkQuaternionRigidTransformGradientDescentOptimizer.h>
#include <itkGradientDescentOptimizer.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkVersorRigid3DTransformOptimizer.h>
#include <itkVersorTransformOptimizer.h>
#include <itkLBFGSOptimizer.h>
#include <itkLBFGSBOptimizer.h>
#include <itkAmoebaOptimizer.h>
#include <itkConjugateGradientOptimizer.h>

#include <itkNormalVariateGenerator.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include <itkCommand.h>


class OptimizerHelper : public itk::Object
{
  
  public:
    typedef OptimizerHelper                             Self;
    typedef itk::Object                                 Superclass;
    typedef itk::SmartPointer< Self >                   Pointer;
    typedef itk::SmartPointer< const Self >             ConstPointer;

    itkTypeMacro( OptimizerHelper, Object );

    itkNewMacro( Self );


    enum OPTIMIZER_TYPE {
      EXHAUSTIVE,
      ONE_PLUS_ONE_EVOLUTIONARY,
      SPSA,
      POWELL,
      POWELL_FRPR,
      GRADIENT_DESCENT,
      GRADIENT_DESCENT_REGULAR_STEP,
      GRADIENT_DESCENT_REGULAR_STEP_VERSOR_RIGID,
      GRADIENT_DESCENT_REGULAR_STEP_VERSOR,
      ROBBINS_MONRO,
      VNL_LBFGS,
      VNL_LBFGSB,
      VNL_AMOEBA,
      VNL_CONJUGATE_GRADIENT
    };


    typedef itk::SingleValuedNonLinearOptimizer                        GenericType;
    typedef itk::ExhaustiveOptimizer                                   ExhaustiveType;
    typedef itk::OnePlusOneEvolutionaryOptimizer                       OnePlusOneType;
    typedef itk::SPSAOptimizer                                         SPSAType;
    typedef itk::PowellOptimizer                                       PowellType;
    typedef itk::FRPROptimizer                                         PowellFRPRType;
    typedef itk::QuaternionRigidTransformGradientDescentOptimizer      QuaternionRigidGDType;
    typedef itk::GradientDescentOptimizer                              GDType;
    typedef itk::RegularStepGradientDescentOptimizer                   RegularStepGDType;
    typedef itk::VersorRigid3DTransformOptimizer                       VersorRigidType;
    typedef itk::VersorTransformOptimizer                              VersorType;
    typedef itk::RobbinsMonroOptimizer                                 RMType;
    typedef itk::LBFGSOptimizer                                        LBFGSType;
    typedef itk::LBFGSBOptimizer                                       LBFGSBType;
    typedef itk::AmoebaOptimizer                                       AmoebaType;
    typedef itk::ConjugateGradientOptimizer                            ConjugateGradientType;
    
    typedef std::vector< boost::any >                                  ParametersType;
    typedef itk::Array< double >                                       TransformParametersType;
    typedef GenericType::ScalesType                                    ScalesType;
    
    void SetType( OPTIMIZER_TYPE type );
    void SetScales( ScalesType scales );
    
    unsigned char GetNumberOfParameters();
    
    void SetParameters( ParametersType param );
    
    void SetOutputStream( std::ostream* os ) { m_os = os; }
    
    ParametersType GetParameters()
    { return m_Parameters; }

    
    itkGetConstMacro(EnableDefaultObserver, bool);
    itkSetMacro(EnableDefaultObserver, bool);
    
    GenericType::Pointer GetOptimizer()
    { return m_Optimizer; }
    
    template <class TRegistration>
    void Connect( TRegistration* r);
    
    template <class T>
    void SetIterationObserverFunction( T* object, typename  itk::MemberCommand<T>::TMemberFunctionPointer memberFunction );
    
    void SPSAGuessParameters();
    
  protected:
    OptimizerHelper();
    ~OptimizerHelper() {};
    
    void DefaultObserver();

    /**
     * Print contents of an OptimizerHelper
     */
    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;
    
  private:
    OptimizerHelper(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    typedef itk::SimpleMemberCommand< Self >  ObserverType;
    
    void CreateOptimizer();

    GenericType::Pointer                m_Optimizer;
    OPTIMIZER_TYPE                      m_Type;
    ParametersType                      m_Parameters;
    ScalesType                          m_Scales;
    bool                                m_EnableDefaultObserver;
    ObserverType::Pointer               m_Observer;
    
    
    ExhaustiveType::Pointer             m_Exhaustive;
    OnePlusOneType::Pointer             m_1plus1;
    SPSAType::Pointer                   m_Spsa;
    PowellType::Pointer                 m_Powell;
    PowellFRPRType::Pointer             m_PowellFRPR;
    QuaternionRigidGDType::Pointer      m_QuatRigid;
    GDType::Pointer                     m_GDescent;
    RegularStepGDType::Pointer          m_RegStepGD;
    VersorRigidType::Pointer            m_VersorRigid;
    VersorType::Pointer                 m_Versor;
    RMType::Pointer                     m_RM;
    LBFGSType::Pointer                  m_LBFGS;
    LBFGSBType::Pointer                 m_LBFGSB;
    AmoebaType::Pointer                 m_Amoeba;
    ConjugateGradientType::Pointer      m_ConjGradient;

    std::ostream*                       m_os;
};

#include "OptimizerHelper.txx"


#endif // OPTIMIZERHELPER_H
