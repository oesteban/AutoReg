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

#ifndef OPTIMIZERHELPER_TXX
#define OPTIMIZERHELPER_TXX

template <class T>
void OptimizerHelper::SetIterationObserverFunction( T* object,  
             typename itk::MemberCommand<T>::TMemberFunctionPointer memberFunction )
{
  ObserverType::Pointer o = ObserverType::New();
  o->SetCallbackFunction( object, & memberFunction );
  o->AddObserver( itk::IterationEvent(), o );
}


template <class TRegistration>
void OptimizerHelper::Connect( TRegistration* r){
    switch( m_Type )
  {
    case EXHAUSTIVE:
      r->SetOptimizer( m_Exhaustive );
      break;
    case ONE_PLUS_ONE_EVOLUTIONARY:
      r->SetOptimizer( m_1plus1 );
      break;
    case SPSA:
      r->SetOptimizer( m_Spsa );
      break;
    case POWELL:
      r->SetOptimizer( m_Powell );
      break;
    case POWELL_FRPR:
      r->SetOptimizer( m_PowellFRPR );
      break;
    case GRADIENT_DESCENT_REGULAR_STEP:
      r->SetOptimizer( m_RegStepGD );
      break;
    case GRADIENT_DESCENT_REGULAR_STEP_VERSOR_RIGID:
      r->SetOptimizer( m_VersorRigid );
      break;
    case GRADIENT_DESCENT_REGULAR_STEP_VERSOR:
      r->SetOptimizer( m_Versor );
      break;
    case ROBBINS_MONRO:
      r->SetOptimizer( m_RM );
      break;
    case VNL_LBFGS:
      r->SetOptimizer( m_LBFGS );
      break;
    case VNL_LBFGSB:
      r->SetOptimizer( m_LBFGSB );
      break;
    case VNL_AMOEBA:
      r->SetOptimizer( m_Amoeba );
      break;
    case VNL_CONJUGATE_GRADIENT:
      r->SetOptimizer( m_ConjGradient );
      break;
    case GRADIENT_DESCENT:
       r->SetOptimizer( m_GDescent );
      break;
    default:
      r->SetOptimizer( m_Optimizer );
      break;
  }
  
  
//   r->SetOptimizer( m_Optimizer );
};


#endif // OPTIMIZERHELPER_TXX
