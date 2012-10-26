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

#include "OptimizerHelper.h"
#include <exception>

OptimizerHelper::OptimizerHelper() :
m_EnableDefaultObserver(true), m_os(&std::cout) {
	m_Observer = ObserverType::New();
	m_Observer->SetCallbackFunction(this, &Self::DefaultObserver);
}

void OptimizerHelper::CreateOptimizer() {
	double maxStep = 1;
	m_Parameters.resize(0);
	switch (m_Type) {
	case EXHAUSTIVE:
		m_Parameters.resize(2);
		m_Parameters[1] = 1.0;
		m_Exhaustive = ExhaustiveType::New();
		m_Optimizer = m_Exhaustive;
		break;
	case ONE_PLUS_ONE_EVOLUTIONARY:
		m_Parameters.resize(6);
		m_1plus1 = OnePlusOneType::New();
		m_Optimizer = m_1plus1;
		break;
	case SPSA:
		m_Parameters.resize(14);
		m_Spsa = SPSAType::New();
		//       m_Spsa->DebugOn();
		m_Optimizer = m_Spsa;
		// Set default parameters: #maxIterations, MaximizeOn, a, A, c, #minIterations, Alpha, Gamma, #perturbations, convDecayRate, Tolerance, Sa, Sc
		m_Parameters[0] = 500u;
		m_Parameters[1] = true;
		m_Parameters[2] = 10000.0;
		m_Parameters[3] = 500 * 0.10;
		m_Parameters[4] = 40.0;
		m_Parameters[9] = 0.0;
		//m_Parameters[10]= 15e-02;
		break;
	case POWELL:
		m_Powell = PowellType::New();
		m_Optimizer = m_Powell;
		break;
	case POWELL_FRPR:
		m_PowellFRPR = PowellFRPRType::New();
		m_Optimizer = m_PowellFRPR;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP:
		m_Parameters.resize(5);
		m_RegStepGD = RegularStepGDType::New();
		m_Optimizer = m_RegStepGD;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR_RIGID:
		m_Parameters.resize(5);
		m_VersorRigid = VersorRigidType::New();
		m_Optimizer = m_VersorRigid;

		// Set default parameters:
		m_Parameters[0] = 1000u;
		m_Parameters[1] = true;
		m_Parameters[2] = 0.9;
		m_Parameters[3] = maxStep;
		m_Parameters[4] = 0.001 * maxStep;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR:
		m_Versor = VersorType::New();
		m_Optimizer = m_Versor;
		break;
	case VNL_LBFGS:
		m_LBFGS = LBFGSType::New();
		m_Optimizer = m_LBFGS;
		break;
	case VNL_LBFGSB:
		m_LBFGSB = LBFGSBType::New();
		m_Optimizer = m_LBFGSB;
		break;
	case VNL_AMOEBA:
		m_Amoeba = AmoebaType::New();
		m_Optimizer = m_Amoeba;
		break;
	case VNL_CONJUGATE_GRADIENT:
		m_ConjGradient = ConjugateGradientType::New();
		m_Optimizer = m_ConjGradient;
		break;
	case ROBBINS_MONRO:
		m_Parameters.resize(5);
		m_RM = RMType::New();
		m_Optimizer = m_RM;
		break;
	default:
		m_Parameters.resize(3);
		m_GDescent = GDType::New();
		m_Optimizer = m_GDescent;
		m_Parameters[0] = 300u;
		m_Parameters[1] = true;
		m_Parameters[2] = 15.0;
		break;
	}

	SetParameters(m_Parameters);

	if (m_EnableDefaultObserver)
		m_Optimizer->AddObserver(itk::IterationEvent(), m_Observer);
}
;

void OptimizerHelper::SetType(OptimizerHelper::OPTIMIZER_TYPE type) {
	m_Type = type;
	CreateOptimizer();
}
;
void OptimizerHelper::SetScales(OptimizerHelper::ScalesType scales) {
	m_Scales = scales;
	m_Optimizer->SetScales(scales);
}
;

unsigned char OptimizerHelper::GetNumberOfParameters() {
	return m_Parameters.size();
}
;

void OptimizerHelper::SetParameters(OptimizerHelper::ParametersType param) {
	bool isMersenneGenerator = false;

	switch (m_Type) {
	case EXHAUSTIVE:
		try {
			m_Exhaustive->SetNumberOfSteps(boost::any_cast<
					ExhaustiveType::StepsType>(param.at(0)));
		} catch (...) {
		}
		try {
			m_Exhaustive->SetStepLength(boost::any_cast<double>(param.at(1)));
		} catch (...) {
		}
		break;
	case ONE_PLUS_ONE_EVOLUTIONARY:
		// MaximumIteration, MaximizeOn, NormalVariateGenerator, Epsilon, InitialRadius, GrowthFactor
		m_1plus1->SetMaximumIteration(
				boost::any_cast<unsigned int>(param.at(0)));
		m_1plus1->SetMaximize(boost::any_cast<bool>(param.at(1)));
		try {
			isMersenneGenerator = boost::any_cast<std::string>(param.at(1))
							== "mersenne";
		} catch (...) {
		}

		if (isMersenneGenerator)
			m_1plus1->SetNormalVariateGenerator(
					itk::Statistics::MersenneTwisterRandomVariateGenerator::New());
		else
			m_1plus1->SetNormalVariateGenerator(
					itk::Statistics::NormalVariateGenerator::New());

		try {
			m_1plus1->SetEpsilon(boost::any_cast<double>(param.at(2)));
		} catch (...) {
		}

		try {
			m_1plus1->SetInitialRadius(boost::any_cast<double>(param.at(3)));
		} catch (...) {
		}

		try {
			m_1plus1->SetGrowthFactor(boost::any_cast<double>(param.at(4)));
		} catch (...) {
		}

		break;
	case SPSA:
		// #maxIterations, MaximizeOn, a, A, c, #minIterations, Alpha, Gamma, #perturbations, convDecayRate, Tolerance, Sa, Sc
		m_Spsa->SetMaximumNumberOfIterations(boost::any_cast<unsigned int>(
				param.at(0)));
		m_Spsa->SetMaximize(boost::any_cast<bool>(param.at(1)));
		try {
			m_Spsa->Seta(boost::any_cast<double>(param.at(2)));
		} catch (...) {
			std::cout << "Warning: SPSA a Optimizer parameter a was not set"
					<< std::endl;
		}
		try {
			m_Spsa->SetA(boost::any_cast<double>(param.at(3)));
		} catch (...) {
		}
		try {
			m_Spsa->Setc(boost::any_cast<double>(param.at(4)));
		} catch (...) {
		}
		try {
			m_Spsa->SetMinimumNumberOfIterations(boost::any_cast<unsigned int>(
					param.at(5)));
		} catch (...) {
		}
		try {
			m_Spsa->SetAlpha(boost::any_cast<double>(param.at(6)));
		} catch (...) {
		}
		try {
			m_Spsa->SetGamma(boost::any_cast<double>(param.at(7)));
		} catch (...) {
		}
		try {
			m_Spsa->SetNumberOfPerturbations(boost::any_cast<unsigned int>(
					param.at(8)));
		} catch (...) {
		}
		try {
			m_Spsa->SetStateOfConvergenceDecayRate(boost::any_cast<double>(
					param.at(13)));
		} catch (...) {
		}
		try {
			m_Spsa->SetTolerance(boost::any_cast<double>(param.at(10)));
		} catch (...) {
		}
		try {
			m_Spsa->SetSa(boost::any_cast<double>(param.at(11)));
		} catch (...) {
		}
		try {
			m_Spsa->SetSc(boost::any_cast<double>(param.at(12)));
		} catch (...) {
		}

		try {
			double step = boost::any_cast<double>(param.at(9));
			if (step > 0.0 ){
				ObserverType::Pointer initObserver = ObserverType::New();
				initObserver->SetCallbackFunction(this, &Self::SPSAGuessParameters);
				m_Spsa->AddObserver(itk::StartEvent(), initObserver);
			}
		} catch (std::exception& e) {
			std::cout << "Parameters couldn't be passed:" << std::endl;
			std::cout << e.what() << std::endl;
		}

		break;
	case POWELL:

	case POWELL_FRPR:

	case GRADIENT_DESCENT_REGULAR_STEP:

	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR_RIGID:
		// NumberOfIterations, MaximizeOn, RelaxationFactor, MaximumStepLength, MinimumStepLength
		m_VersorRigid->SetNumberOfIterations(boost::any_cast<unsigned int>(
				param.at(0)));
		m_VersorRigid->SetMaximize(boost::any_cast<bool>(param.at(1)));
		m_VersorRigid->SetRelaxationFactor(boost::any_cast<double>(param.at(2)));
		m_VersorRigid->SetMaximumStepLength(
				boost::any_cast<double>(param.at(3)));
		m_VersorRigid->SetMinimumStepLength(
				boost::any_cast<double>(param.at(4)));
		break;
	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR:

	case VNL_LBFGS:

	case VNL_LBFGSB:

	case VNL_AMOEBA:

	case VNL_CONJUGATE_GRADIENT:

	case ROBBINS_MONRO:
		try {
			m_RM->SetNumberOfIterations(boost::any_cast<unsigned int>(param.at(
					0)));
		} catch (...) {
		}
		try {
			m_RM->SetMaximize(boost::any_cast<bool>(param.at(1)));
		} catch (...) {
		}
		try {
			m_RM->SetA(boost::any_cast<double>(param.at(2)));
		} catch (...) {
		}
		try {
			m_RM->Seta(boost::any_cast<double>(param.at(3)));
		} catch (...) {
		}
		try {
			m_RM->SetAlpha(boost::any_cast<double>(param.at(4)));
		} catch (...) {
		}

		break;
	default:
		// GRADIENT_DESCENT: NumberOfIterations, MaximizeOn, LearningRate
		m_GDescent->SetNumberOfIterations(boost::any_cast<unsigned int>(
				param.at(0)));
		m_GDescent->SetMaximize(boost::any_cast<bool>(param.at(1)));
		m_GDescent->SetLearningRate(boost::any_cast<double>(param.at(2)));
		break;
	}
	m_Parameters=param;
}
;

void OptimizerHelper::SPSAGuessParameters() {
	try {
		m_Spsa->GuessParameters(4u, boost::any_cast<double>(m_Parameters.at(9)));
	} catch (itk::ExceptionObject & err) {
		std::cout << "Parameters couldn't be guessed:" << std::endl;
		std::cout << err << std::endl;
	}
}
;

void OptimizerHelper::DefaultObserver() {
	itk::Optimizer::ParametersType inc;
	static itk::Optimizer::ParametersType lastParameters;
	switch (m_Type) {
	case EXHAUSTIVE:
		inc = m_Exhaustive->GetCurrentPosition() - m_Exhaustive->GetInitialPosition();

		*m_os << std::setw(10) << std::setfill(' ') << std::left << sqrt(inc[0]
		                                                                     * inc[0] + inc[1] * inc[1] + inc[2] * inc[2]) << " ";
		*m_os << std::setw(10) << std::setfill(' ') << std::left << sqrt(inc[3]
		                                                                     * inc[3] + inc[4] * inc[4] + inc[5] * inc[5]) << " ";
		*m_os << std::setw(8) << std::setfill(' ') << std::left
				<< m_Exhaustive->GetCurrentValue();
		*m_os << "[" << inc << "]" << std::endl;
		break;
	case ONE_PLUS_ONE_EVOLUTIONARY:
		*m_os << m_1plus1->GetCurrentIteration() << "   ";
		*m_os << m_1plus1->GetValue() << "   ";
		*m_os << m_1plus1->GetCurrentPosition() << std::endl;
		break;
	case SPSA:
		*m_os << m_Spsa->GetCurrentIteration() << "   ";
		*m_os << std::setw(6) << std::left << m_Spsa->GetValue() << "   ";
		if( lastParameters.GetNumberOfElements() == 0 ) lastParameters= m_Spsa->GetInitialPosition();
		*m_os << "[" << m_Spsa->GetCurrentPosition() - lastParameters << "]";
		lastParameters = m_Spsa->GetCurrentPosition();

		*m_os << "  " << m_Spsa->GetStateOfConvergence() << std::endl;

		/*

		for (unsigned int i = 0;i<m_Spsa->GetGradient().GetSize();i++)
		  *m_os << (m_Spsa->GetGradient()[i] * m_Spsa->GetLearningRate() ) << "  ";
		*m_os << "]";
		*m_os << m_Spsa->GetLearningRate() << "  ";
		*m_os << m_Spsa->GetGradientMagnitude() << "  ";
		*m_os << m_Spsa->GetGradient() << "  ";

		*/
		break;
	case POWELL:
		*m_os << m_Powell->GetCurrentIteration() << "   ";
		*m_os << m_Powell->GetValue() << "   ";
		*m_os << m_Powell->GetCurrentPosition() << std::endl;
		break;
	case POWELL_FRPR:
		*m_os << m_PowellFRPR->GetCurrentIteration() << "   ";
		*m_os << m_PowellFRPR->GetValue() << "   ";
		*m_os << m_PowellFRPR->GetCurrentPosition() << std::endl;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP:
		*m_os << m_RegStepGD->GetCurrentIteration() << "   ";
		*m_os << m_RegStepGD->GetValue() << "   ";
		*m_os << m_RegStepGD->GetCurrentPosition() << std::endl;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR_RIGID:
		*m_os << m_VersorRigid->GetCurrentIteration() << "   ";
		*m_os << m_VersorRigid->GetValue() << "   ";
		*m_os << m_VersorRigid->GetCurrentPosition() << std::endl;
		break;
	case GRADIENT_DESCENT_REGULAR_STEP_VERSOR:
		*m_os << m_Versor->GetCurrentIteration() << "   ";
		*m_os << m_Versor->GetValue() << "   ";
		*m_os << m_Versor->GetCurrentPosition() << std::endl;
		break;
	case VNL_LBFGS:
		//*m_os << m_LBFGS->Get << "   ";
		*m_os << m_LBFGS->GetValue() << "   ";
		*m_os << m_LBFGS->GetCurrentPosition() << std::endl;
		break;
	case VNL_LBFGSB:
		*m_os << m_LBFGSB->GetCurrentIteration() << "   ";
		*m_os << m_LBFGSB->GetValue() << "   ";
		*m_os << m_LBFGSB->GetCurrentPosition() << std::endl;
		break;
	case VNL_AMOEBA:
		*m_os << m_Amoeba->GetValue() << "   ";
		*m_os << m_Amoeba->GetCurrentPosition() << std::endl;
		break;
	case VNL_CONJUGATE_GRADIENT:
		*m_os << m_ConjGradient->GetCurrentIteration() << "   ";
		*m_os << m_ConjGradient->GetValue() << "   ";
		*m_os << m_ConjGradient->GetCurrentPosition() << std::endl;
		break;
	case ROBBINS_MONRO:
		*m_os << m_RM->GetCurrentIteration() << "   ";
		*m_os << std::setw(6) << std::left << m_RM->GetValue() << "   ";
//		*m_os << m_RM->GetCurrentPosition() << "    ";

		if( lastParameters.GetNumberOfElements() == 0 ) lastParameters= m_RM->GetInitialPosition();
		*m_os << "[" << m_RM->GetCurrentPosition() - lastParameters << "]";
		lastParameters = m_RM->GetCurrentPosition();

		*m_os << m_RM->GetLearningRate() << std::endl;
		break;
	default:
		*m_os << m_GDescent->GetCurrentIteration() << "   ";
		*m_os << m_GDescent->GetValue() << "   ";
		*m_os << m_GDescent->GetCurrentPosition() << std::endl;
		break;
	}
}
;

void OptimizerHelper::PrintSelf(std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	m_Optimizer->Print(os, indent);

	os << indent << indent << "* Scales: " << m_Optimizer->GetScales()
					<< std::endl;
}

/*
 int iteration = optimizer->GetCurrentIteration();
 double actualLearningRate = optimizer->GetLearningRate();
 static double bestMetric = 0;
 static int last_index = 0;
 const double metric_value = optimizer->GetValue();

 output << "[" << iteration << "] ";
 output << metric_value << ":: ";
 output << optimizer->GetCurrentPosition()  << "." << std::endl;


 if (last_index > -1 && last_index < OPTIMIZER_ITERATIONS_UPDATE)
 {
 if ( (optimizer->GetMinimize() && metric_value >= bestMetric) ||
 (!optimizer->GetMinimize() && metric_value <= bestMetric) )
 {
 last_index++;
 if (last_index == OPTIMIZER_ITERATIONS_UPDATE)
 {
 double newLR = actualLearningRate*(0.60 * (1.0 - (iteration/optimizer->GetNumberOfIterations())) );
 if ( newLR > GRADIENT_DESCENT_MIN_LR )
 {
 optimizer->SetLearningRate( newLR );
 }
 else
 {
 optimizer->StopOptimization();
 bestMetric = 0;
 }
 last_index= 0;
 }
 }
 else
 {
 last_index = 0;
 bestMetric = metric_value;
 }
 }*/
