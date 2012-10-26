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
#ifndef TRANSFORMHELPER_TXX
#define TRANSFORMHELPER_TXX

#include "TransformHelper.h"

#include <itkTransformFileWriter.h>
#include <itkImageRegionConstIterator.h>

#include <itkShrinkImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageDuplicator.h>

template<class TFixedImage, class TMovingImage, class TScalarType>
TransformHelper<TFixedImage, TMovingImage, TScalarType>::TransformHelper() :
	m_InitType(INIT_NONE), m_SwapLevel(false) {
	
	
}


template<class TFixedImage, class TMovingImage, class TScalarType>
void TransformHelper<TFixedImage, TMovingImage, TScalarType>::SetTransformType( typename TransformHelper<TFixedImage, TMovingImage, TScalarType>::TRANSFORM_TYPE type ) { 
	m_Type = type;
	
	switch (m_Type) {
		case EULER_3D:
			m_EulerTransform = Instantiate<EulerType> ();
			break;
		case EULER_3D_CENTERED:
			m_CenteredEulerTransform = Instantiate<CenteredEulerType> ();
			break;
		case VERSOR_RIGID_3D:
			m_VersorRigidTransform = Instantiate<VersorRigidType> ();
			break;
		case VERSOR_NON_RIGID_3D:
			m_VersorTransform = Instantiate<VersorType> ();
			break;
		case AFFINE_7P:
			m_Affine7PTransform = Instantiate<Affine7PType> ();
			break;
		case AFFINE_9P:
			m_Affine9PTransform = Instantiate<Affine9PType> ();
			break;
		case AFFINE_15P:
			m_Affine15PTransform = Instantiate<Affine15PType> ();
			break;
		case BSPLINE:
			std::cout << "BSpline transform Not Implemented" << std::endl;
			break;
		default:
			// TODO throw exception
			std::cout << "Unsoported transform code" << std::endl;
			break;
	}
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void TransformHelper<TFixedImage, TMovingImage, TScalarType>::PrintSelf(
		std::ostream& os, itk::Indent indent) const {
	Superclass::PrintSelf(os, indent);

	if (m_InitTransform.IsNotNull()) {
		os << indent << " * m_InitTransform: " << std::endl;
		m_InitTransform->Print(os, indent);
	}

	if (m_Transform.IsNotNull()) {
		os << indent << " * m_Transform: " << std::endl;
		m_Transform->Print(os, indent);
	}
}

template<class TFixedImage, class TMovingImage, class TScalarType>
typename TransformHelper<TFixedImage, TMovingImage, TScalarType>::GenericType::Pointer TransformHelper<
		TFixedImage, TMovingImage, TScalarType>::GetInverseTransform() {
	typedef typename RigidType::Pointer TFormPointer;

	TFormPointer inverse = RigidType::New();

	TFormPointer rigid_tf = dynamic_cast<RigidType*> (m_Transform.GetPointer());
	inverse->SetCenter(rigid_tf->GetCenter());
	rigid_tf->GetInverse(inverse);

	typename GenericType::Pointer result(inverse.GetPointer());

	return result;
}

template<class TFixedImage, class TMovingImage, class TScalarType>
void TransformHelper<TFixedImage, TMovingImage, TScalarType>::SaveTransformToFile(
		std::string fileName) {
	typedef itk::TransformFileWriter TFWriter;
	TFWriter::Pointer tw = TFWriter::New();
	tw->SetFileName(fileName);
	tw->AddTransform(m_Transform);
	tw->SetAppendOn();
	tw->Update();
}

template<class TFixedImage, class TMovingImage, class TScalarType>
template<class TTransform>
typename TTransform::Pointer TransformHelper<TFixedImage, TMovingImage,
		TScalarType>::Instantiate() {
	typename TTransform::Pointer transform = TTransform::New();
	transform->SetIdentity();
	m_NumberOfParameters = transform->GetNumberOfParameters();
	//m_Transform = transform;
	return transform;
}

template<class TFixedImage, class TMovingImage, class TScalarType>
template<class TTransform>
typename TTransform::Pointer TransformHelper<TFixedImage, TMovingImage,
		TScalarType>::InitializeTransform() {
	/*
	typename TTransform::Pointer transform = TTransform::New();
	transform->SetIdentity();
	m_Transform = transform;
	*/
	typename TTransform::Pointer transform = this->Instantiate< TTransform > ();
	m_Transform = transform;

	if ((m_InitTransform.IsNull() && (m_FixedImage.IsNull()
			|| m_MovingImage.IsNull())) || (m_InitType == INIT_NONE)) {
		m_InitType = INIT_NONE;
		return transform;
	}


	// Transform initialization. If init transform is present, set it as initialization AND finish.
	if (m_InitTransform.IsNotNull()) {
		// If transform is equivalent, copy it
		if (transform->GetParameters().size() == m_InitTransform->GetParameters().size())
			transform = static_cast<TTransform *> (m_InitTransform.GetPointer());

		// If not, set manually the parameters
		else {
			typename GenericType::ParametersType init_p = m_InitTransform->GetParameters();
			// unsigned int numOfParameters = init_p.size();

			transform->SetFixedParameters(m_InitTransform->GetFixedParameters());
			typename TTransform::MatrixType m;
			m.SetIdentity();

			for (unsigned int i = 0; i < 3; i++)
				for (unsigned int j = 0; j < 3; j++) m(i, j) = init_p[i * 3 + j];

			transform->SetMatrix(m);

			typename TTransform::OutputVectorType t;
			for (unsigned int k = 0; k < 3; k++) t[k] = init_p[k + 9];

			transform->SetTranslation(t);
		}
		m_Transform = transform;
		return transform;
	}


	// Other initializations
	if (m_InitType == INIT_AFFINE || m_InitType == INIT_PRINCIPAL_AXIS ) {
		typename TTransform::ParametersType param( transform->GetNumberOfParameters() );
		param.Fill(0.0);

		typedef itk::ImageMomentsCalculator<FixedImageType> FixedMomentsCalc;
		typename FixedMomentsCalc::Pointer fixedCalculator =
				FixedMomentsCalc::New();
		fixedCalculator->SetImage(m_FixedImage);
		fixedCalculator->Compute();

		typedef itk::ImageMomentsCalculator<MovingImageType> MovingMomentsCalc;
		typename MovingMomentsCalc::Pointer movingCalculator =
				MovingMomentsCalc::New();
		movingCalculator->SetImage(m_MovingImage);
		movingCalculator->Compute();

		typename FixedMomentsCalc::VectorType fixedCenter =
				fixedCalculator->GetCenterOfGravity();
		typename MovingMomentsCalc::VectorType movingCenter =
				movingCalculator->GetCenterOfGravity();

		typename TTransform::OutputVectorType translationVector;
		typename TTransform::InputPointType rotationCenter;

		for (unsigned int i = 0; i < FixedImageDimension; i++) {
			rotationCenter[i] = fixedCenter[i];
			translationVector[i] = movingCenter[i] - fixedCenter[i];
			param[i + 3] = translationVector[i];
		}

		transform->SetCenter(rotationCenter);
		transform->SetTranslation( translationVector );

		if ( m_InitType == INIT_PRINCIPAL_AXIS ) {
			//typename FixedMomentsCalc::AffineTransformPointer f_pa = fixedCalculator->GetPrincipalAxesToPhysicalAxesTransform();
			//typename FixedMomentsCalc::AffineTransformType::MatrixType m_fixed = tf_pa->GetMatrix();
/*			typename FixedMomentsCalc::AffineTransformType::MatrixType m_fixed = fixedCalculator->GetCentralMoments();

			typename FixedMomentsCalc::AffineTransformType::InputVectorType v_fixed[3];
			for (unsigned int i = 0; i< FixedImageDimension; i++) {
				for (unsigned int j = 0; j< FixedImageDimension; j++) {
					v_fixed[i][j] = m_fixed(i,j);
				}
				std::cout << "Fixed Vector [" << i << "]=" << v_fixed[i] << ", Norm=" << v_fixed[i].GetNorm() << std::endl;
				v_fixed[i].Normalize();
			}


			//typename MovingMomentsCalc::AffineTransformPointer m_pa = movingCalculator->GetPrincipalAxesToPhysicalAxesTransform();
			//typename MovingMomentsCalc::AffineTransformType::MatrixType m_moving = pa_tf->GetMatrix();
			typename MovingMomentsCalc::AffineTransformType::MatrixType m_moving = movingCalculator->GetCentralMoments();

			typename MovingMomentsCalc::AffineTransformType::InputVectorType v_moving[3];
			for (unsigned int i = 0; i< MovingImageDimension; i++) {
				for (unsigned int j = 0; j< MovingImageDimension; j++) {
					v_moving[i][j] = m_moving(i,j);
				}
				std::cout << "Moving Vector [" << i << "]=" << v_moving[i] << ", Norm=" << v_moving[i].GetNorm()<< std::endl;
				v_moving[i].Normalize();
			}

			typedef itk::Euler3DTransform< double > TGenerator;
			TGenerator::Pointer t = TGenerator::New();
			t->SetIdentity();
*/
/*
			for ( unsigned int i = 0; i<3; i++) {
				for ( unsigned int axis = 0; axis<3; axis++) {
					TGenerator::InputVectorType axisV; axisV.Fill(0.0); axisV[axis]=1.0;
					TGenerator::InputVectorType axisF; axisF.Fill(0.0); axisF[axis]=v_fixed[i][axis];
					TGenerator::InputVectorType axisM; axisM.Fill(0.0); axisM[axis]=v_moving[i][axis];
					std::cout << "Angle fixed[" << i<< "]" << "- axis[" << axis << "]=" << 	acos( (double) (axisF * axisV) ) << std::endl;
					std::cout << "Angle moving[" << i<< "]" << "- axis[" << axis << "]=" << 	acos( (double) (axisM * axisV) ) << std::endl;

					std::cout << "Angle fixed-moving[" << i<< "] - axis[="<< axis << "]=" <<  	acos( (double) (axisF * axisM) ) << std::endl;
				}
			}*/
/*
			TGenerator::InputVectorType axisV; axisV.Fill(0.0); axisV[2]=1.0;
			TGenerator::InputVectorType axisF; axisF.Fill(0.0); axisF[2]=v_fixed[2][2];
			TGenerator::InputVectorType axisM; axisM.Fill(0.0); axisM[2]=v_moving[2][2];

			t->SetRotation( acos( (double) (axisF * axisV) ) + acos( (double) (axisM * axisV) ), 0.0, 0.0 );
*/

			//transform->SetMatrix( t->GetMatrix() );
			transform->SetMatrix(  movingCalculator->GetPrincipalAxes()* (fixedCalculator->GetPrincipalAxes()).GetInverse() );
			//transform->SetMatrix(  (fixedCalculator->GetPrincipalAxes()).GetInverse() );

		}
		else if ( m_InitType == INIT_AFFINE ) { 			// Isometric Scaling depending on Mass
			if (param.GetSize() > 6) {
				double fixedMass = fixedCalculator->GetTotalMass();
				double movingMass = movingCalculator->GetTotalMass();
				double fixedVol = 1.0;
				double movingVol = 1.0;

				for (unsigned i = 0; i < FixedImageDimension; i++) {
					fixedVol *= m_FixedImage->GetSpacing()[i];
				}
				for (unsigned i = 0; i < MovingImageDimension; i++) {
					movingVol *= m_MovingImage->GetSpacing()[i];
				}

				double scale = std::pow((movingMass * movingVol) / (fixedMass
						* fixedVol), 1.0 / 3.0);
				for (unsigned i = 6; i < param.GetSize(); i++) {
					param[i] = scale;
				}
			}
			transform->SetParameters(param);
		}
	}

	typedef itk::CenteredTransformInitializer<TTransform, FixedImageType,
			MovingImageType> Initializer;
	typename Initializer::Pointer init;

	if (m_InitType == INIT_GEOMETRY || m_InitType == INIT_MOMENTS) {
		init = Initializer::New();
		init->SetTransform(transform);
		init->SetFixedImage(m_FixedImage);
		init->SetMovingImage(m_MovingImage);

		if (m_InitType == INIT_GEOMETRY)
			init->GeometryOn();
		else
			init->MomentsOn();
	}

	if (init.IsNotNull())
		init->InitializeTransform();

	if (m_SwapLevel)
		transform->GetInverse(transform);

	return transform;
}
;

template<class TFixedImage, class TMovingImage, class TScalarType>
template<class TInputImage, class TReferenceImage>
TReferenceImage* TransformHelper<TFixedImage, TMovingImage, TScalarType>::ResampleImage(
		typename TInputImage::ConstPointer input,
		typename TReferenceImage::ConstPointer ref, bool useInverse,
		std::string fileName) {
	typedef itk::ResampleImageFilter<TInputImage, TReferenceImage> Resampler;
	typename Resampler::Pointer resample = Resampler::New();
	resample->SetInput(input);
	resample->SetUseReferenceImage(true);
	resample->SetReferenceImage(ref);
	resample->SetDefaultPixelValue(0);

	if (useInverse) {
		resample->SetTransform(GetInverseTransform());
	} else {
		resample->SetTransform(m_Transform);
	}

	if (fileName.compare("") != 0)
		SaveImageToFile<TReferenceImage> (resample->GetOutput(), fileName);

	return resample->GetOutput();
}

template<class TFixedImage, class TMovingImage, class TScalarType>
template<class TRegistration>
void TransformHelper<TFixedImage, TMovingImage, TScalarType>::Connect(
		TRegistration* r) {
	typename VersorRigidType::ParametersType p;
	switch (m_Type) {
	case EULER_3D:
		m_EulerTransform = InitializeTransform<EulerType> ();
		r->SetTransform(m_EulerTransform);
		r->SetInitialTransformParameters(m_EulerTransform->GetParameters());
		break;
	case EULER_3D_CENTERED:
		m_CenteredEulerTransform = InitializeTransform<CenteredEulerType> ();
		r->SetTransform(m_CenteredEulerTransform);
		r->SetInitialTransformParameters(
				m_CenteredEulerTransform->GetParameters());
		break;
	case VERSOR_RIGID_3D:
		m_VersorRigidTransform = InitializeTransform<VersorRigidType> ();
		p = m_VersorRigidTransform->GetParameters();
		r->SetTransform(m_VersorRigidTransform);
		r->SetInitialTransformParameters(
				m_VersorRigidTransform->GetParameters());
		break;
	case VERSOR_NON_RIGID_3D:
		m_VersorTransform = InitializeTransform<VersorType> ();
		r->SetTransform(m_VersorTransform);
		r->SetInitialTransformParameters(m_VersorTransform->GetParameters());
		break;
	case AFFINE_7P:
		m_Affine7PTransform = InitializeTransform<Affine7PType> ();
		r->SetTransform(m_Affine7PTransform);
		r->SetInitialTransformParameters(m_Affine7PTransform->GetParameters());
		break;
	case AFFINE_9P:
		m_Affine9PTransform = InitializeTransform<Affine9PType> ();
		r->SetTransform(m_Affine9PTransform);
		r->SetInitialTransformParameters(m_Affine9PTransform->GetParameters());
		break;
	case AFFINE_15P:
		m_Affine15PTransform = InitializeTransform<Affine15PType> ();
		r->SetTransform(m_Affine15PTransform);
		r->SetInitialTransformParameters(m_Affine15PTransform->GetParameters());
		break;
	case BSPLINE:
		std::cout << "BSpline transform Not Implemented" << std::endl;
		break;
	default:
		r->SetTransform(m_Transform);
		r->SetInitialTransformParameters(m_Transform->GetParameters());
		break;
	}
}
;

#endif // TRANSFORMHELPER_TXX
