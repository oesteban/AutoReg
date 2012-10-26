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

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkVersorRigid3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkGIBUB3DTransform.h>
#include <itkOrientImageFilter.h>

#include <itkImageMomentsCalculator.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryContourImageFilter.h>
#include <itkLineConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkConstantBoundaryCondition.h>
#include <itkImageDuplicator.h>

#include <itkContinuousIndex.h>
#include <itkImageRegionIterator.h>

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0

inline double median(std::vector<double> &v) {
	int n = v.size() * 0.5;
	std::nth_element(v.begin(), v.begin() + n, v.end());
	return v[n];
}

inline double max(std::vector<double> &v) {
	double maxvalue = *std::max_element(v.begin(), v.end());
	return maxvalue;
}

inline double mean(std::vector<double> &v) {
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	return sum / v.size();
}

int main(int argc, char *argv[]) {
	// Parameters assignment
	unsigned int case_id;
	// Images names
	std::string bas_im_name, act_im_name, ref_im_name, brain_msk_name, les_msk_name;
	// Transforms files
	std::string bas_tfs_file, act_tfs_file, reg_tfs_file, coreg_tfs_file;
	// Other parameters
	std::string outputPrefix;
	bool output_all, output_msk, output_ref, output_bas, output_act, output_gld, output_dummy;
	bool lesionMask = false;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("case-id", bpo::value<unsigned int>(&case_id)->required(), "test case id")
			("II", bpo::value<std::string>(&bas_im_name)->required(),"basal (interictal) image file")
			("I", bpo::value<std::string>(&act_im_name)->required(), "activated (ictal) image file")
			("ref,R", bpo::value<std::string>(&ref_im_name)->required(),"reference image file")
			("brain,M", bpo::value<std::string>(&brain_msk_name), "skull-stripped brain image file")
			("lesion,L", bpo::value<std::string>(&les_msk_name),"skull-stripped brain image file")
			("bas-tfm,b", bpo::value<std::string>(&bas_tfs_file)->required(),"Gold transforms file, from ref to basal space. itk::VersorRigid3DTranform transforms in ITK format.")
			("act-tfm,a", bpo::value<std::string>(&act_tfs_file)->required(),"Gold transforms file, from ref to ictal space. itk::VersorRigid3DTranform transforms in ITK format.")
			("reg-tfm,r", bpo::value<std::string>(&reg_tfs_file),"Registration test transforms file, from interictal to ictal space. itk::VersorRigid3DTranform transforms in ITK format.")
			("coreg-tfm,c",bpo::value<std::string>(&coreg_tfs_file)->required(), "Co-Registration test transforms file, from interictal to ref space. itk::VersorRigid3DTranform transforms in ITK format.")
			("output-prefix,o", bpo::value<std::string>(&outputPrefix)->default_value("./"), "prefix for output files")
			("output-all", bpo::bool_switch(&output_all)->default_value(false),	"output resampled images in all spaces")
			("output-ref", bpo::bool_switch(&output_ref)->default_value(false), "output resampled images in reference space")
			("output-bas", bpo::bool_switch(&output_bas)->default_value(false), "output resampled images in interictal space")
			("output-act", bpo::bool_switch(&output_act)->default_value(false), "output resampled images in ictal")
			("output-gld", bpo::bool_switch(&output_gld)->default_value(false), "output resampled images in ref space using gold transforms")
			("output-dummy",bpo::bool_switch(&output_dummy)->default_value(false), "output dummy resampled images")
			("output-masks", bpo::bool_switch(&output_msk)->default_value(false), "output masks");

	bpo::positional_options_description pdesc;
	pdesc.add("case-id", 1);
	bpo::variables_map vmap;
	bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).style(bpo::command_line_style::default_style | bpo::command_line_style::allow_long_disguise).run(), vmap);

	if (vmap.count("help")) {
		std::cout << desc << std::endl;
		return 1;
	}

	try {
		bpo::notify(vmap);
	} catch (boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::required_option> > &err) {
		std::cout << "Error: " << err.what() << std::endl;
		std::cout << desc << std::endl;
		return EXIT_FAILURE;
	}

	std::stringstream s_id;
	s_id << case_id;

	typedef float PixelRange;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelRange, Dimension> Image;
	typedef itk::ImageRegionIterator<Image> Iterator;
	typedef itk::Euler3DTransform<double> EulerTransform;
	typedef itk::VersorRigid3DTransform<double> Transform;
	typedef itk::ImageFileReader<Image> Reader;
	typedef itk::ImageFileWriter<Image> Writer;
	typedef itk::ResampleImageFilter<Image, Image> Resampler;
	typedef itk::ImageDuplicator<Image> Duplicator;
	typedef itk::ImageMomentsCalculator<Image> MomentsCalculator;

	// ----------------------------------------------------------------
	// LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------
	Reader::Pointer r_bas = Reader::New();
	Reader::Pointer r_act = Reader::New();
	Reader::Pointer r_ref = Reader::New();
	Reader::Pointer r_les = Reader::New();
	Reader::Pointer r_brain = Reader::New();
	Writer::Pointer w = Writer::New();

	r_bas->SetFileName(bas_im_name);
	r_act->SetFileName(act_im_name);
	r_ref->SetFileName(ref_im_name);
	r_les->SetFileName(les_msk_name);
	r_brain->SetFileName(brain_msk_name);

	try {
		r_bas->Update();
		r_act->Update();
		r_ref->Update();
		r_brain->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	try {
		r_les->Update();
		lesionMask = true;
	} catch (...) {
		lesionMask = false;
	}

	// Set images pointers
	Image::Pointer i_bas = r_bas->GetOutput();
	Image::Pointer i_act = r_act->GetOutput();
	Image::Pointer i_ref = r_ref->GetOutput();
	Image::Pointer i_mask = r_les->GetOutput();
	Image::Pointer i_brain = r_brain->GetOutput();
/*
	std::cout << "i_bas -----------------------------------------------------------------------------------------" << std::endl;
	i_bas->Print( std::cout );
	std::cout << "i_act -----------------------------------------------------------------------------------------" << std::endl;
	i_act->Print( std::cout );
	std::cout << "i_ref -----------------------------------------------------------------------------------------" << std::endl;
	i_ref->Print( std::cout );
	std::cout << "i_mask -----------------------------------------------------------------------------------------" << std::endl;
	i_mask->Print( std::cout );
	std::cout << "i_brain -----------------------------------------------------------------------------------------" << std::endl;
	i_brain->Print( std::cout );
*/
	Duplicator::Pointer dup = Duplicator::New();
	dup->SetInputImage(i_brain);
	dup->Update();
	Image::Pointer mask_out = dup->GetOutput();

	// Compute alignment of centers vector
	Image::SizeType s_bas = i_bas->GetLargestPossibleRegion().GetSize();
	Image::PointType o_bas = i_bas->GetOrigin();
	Image::SpacingType sp_bas = i_bas->GetSpacing();
	Image::SizeType s_act = i_act->GetLargestPossibleRegion().GetSize();
	Image::PointType o_act = i_act->GetOrigin();
	Image::SpacingType sp_act = i_act->GetSpacing();
	Image::SizeType s_ref = i_ref->GetLargestPossibleRegion().GetSize();
	Image::PointType o_ref = i_ref->GetOrigin();
	Image::SpacingType sp_ref = i_ref->GetSpacing();

	Image::PointType center_bas;
	double center_index_bas[3];
	center_index_bas[0] = (s_bas[0] - 1) * 0.5;
	center_index_bas[1] = (s_bas[1] - 1) * 0.5;
	center_index_bas[2] = (s_bas[2] - 1) * 0.5;
	i_bas->TransformContinuousIndexToPhysicalPoint<double> (center_index_bas, center_bas);

	Image::PointType center_ref;
	double center_index_ref[3];
	center_index_ref[0] = (s_ref[0] - 1) * 0.5;
	center_index_ref[1] = (s_ref[1] - 1) * 0.5;
	center_index_ref[2] = (s_ref[2] - 1) * 0.5;
	i_ref->TransformContinuousIndexToPhysicalPoint<double> (center_index_ref, center_ref);

	Image::PointType center_act;
	double center_index_act[3];
	center_index_act[0] = (s_act[0] - 1) * 0.5;
	center_index_act[1] = (s_act[1] - 1) * 0.5;
	center_index_act[2] = (s_act[2] - 1) * 0.5;
	i_act->TransformContinuousIndexToPhysicalPoint<double> (center_index_act, center_act);

	Image::PointType o_zero;
	o_zero.Fill(0.0);

	Transform::InputVectorType ov_rm2ii = o_ref - o_bas;
	Transform::InputVectorType ov_rm2i = o_ref - o_act;
	Transform::InputVectorType ov_ii2i = o_bas - o_act;

	Transform::InputVectorType v_rm2ii = center_ref - center_bas;
	Transform::InputVectorType v_rm2i = center_ref - center_act;
	Transform::InputVectorType v_ii2i = center_bas - center_act;

	Image::PointType bas_ref_orig = o_zero - v_rm2ii;

	// ---------------------------------------------------------------------------------------
	// TEST TRANSFORMS CONFIGURATION ---------------------------------------------------------
	// ---------------------------------------------------------------------------------------

	Transform::Pointer tf_reg;
	Transform::Pointer tf_coreg;
	Transform::Pointer tf_reg_i;
	Transform::Pointer tf_coreg_act;

	typedef itk::TransformFileReader::TransformListType * TFList;
	// ---------------------------------------------------------------------------
	// Registration Test transform (Interictal->Ictal): tf_reg_i

	itk::TransformFileReader::Pointer tf_reg_i_reader = itk::TransformFileReader::New();
	tf_reg_i_reader->SetFileName(reg_tfs_file);
	tf_reg_i_reader->Update();

	TFList tfs_reg_list = tf_reg_i_reader->GetTransformList();
	itk::TransformFileReader::TransformListType::const_iterator tfs_reg_it = tfs_reg_list->begin();

	std::advance(tfs_reg_it, (case_id - 1));

	if (!strcmp((*tfs_reg_it)->GetNameOfClass(), "VersorRigid3DTransform")) {
		tf_reg_i = static_cast<Transform *> ((*tfs_reg_it).GetPointer());
	} else {
		std::cout << "Reading registration test transform failed" << std::endl;
	}

	// ---------------------------------------------------------------------------
	// Coregistration Test transform (Interictal->RM): tf_test_bas2rm

	itk::TransformFileReader::Pointer tf_coreg_reader = itk::TransformFileReader::New();
	tf_coreg_reader->SetFileName(coreg_tfs_file);
	tf_coreg_reader->Update();

	TFList tfs_coreg_list = tf_coreg_reader->GetTransformList();
	itk::TransformFileReader::TransformListType::const_iterator tfs_coreg_it = tfs_coreg_list->begin();

	std::advance(tfs_coreg_it, (case_id - 1));

	if (!strcmp((*tfs_coreg_it)->GetNameOfClass(), "VersorRigid3DTransform")) {
		tf_coreg = static_cast<Transform *> ((*tfs_coreg_it).GetPointer());
	} else {
		std::cout << "Reading coregistration test transform failed" << std::endl;
	}

	tf_reg = Transform::New();
	tf_reg->SetIdentity();
	tf_reg->SetCenter(tf_reg_i->GetCenter());
	tf_reg_i->GetInverse(tf_reg);

	// Old stuff ---------------------------------------------------------------------------
	Transform::Pointer tf_test_bas2act = Transform::New();
	tf_test_bas2act->SetParameters(tf_reg_i->GetParameters());
	tf_test_bas2act->SetFixedParameters(tf_reg_i->GetFixedParameters());

	Image::PointType new_origin;
	new_origin.Fill(0.0);
	new_origin[1] = (i_act->GetLargestPossibleRegion().GetSize()[1] - 1) * i_act->GetSpacing()[1];
	Transform::InputVectorType v = i_act->GetOrigin() - new_origin;
	tf_test_bas2act->SetCenter(tf_test_bas2act->GetCenter() + v);
	Transform::Pointer tf_test_bas2act_inv = Transform::New();
	tf_test_bas2act->GetInverse(tf_test_bas2act_inv);

	Transform::Pointer tf_test_bas2rm = Transform::New();
	tf_test_bas2rm->SetParameters(tf_coreg->GetParameters());
	tf_test_bas2rm->SetFixedParameters(tf_coreg->GetFixedParameters());

	Transform::Pointer tf_test_bas2rm_inv = Transform::New();
	tf_test_bas2rm_inv->SetCenter(tf_test_bas2rm->GetCenter());
	tf_test_bas2rm->GetInverse(tf_test_bas2rm_inv);

	Transform::Pointer tf_test_act2rm = Transform::New();
	tf_test_act2rm->SetFixedParameters(tf_test_bas2rm->GetFixedParameters());
	tf_test_act2rm->SetParameters(tf_test_bas2rm->GetParameters());
	tf_test_act2rm->Compose(tf_test_bas2act_inv, true);

	Transform::Pointer tf_test_act2rm_inv = Transform::New();
	tf_test_act2rm_inv->SetCenter(tf_test_act2rm->GetCenter());
	tf_test_act2rm->GetInverse(tf_test_act2rm_inv);

	// New stuff ---------------------------------------------------------------------------
	Transform::MatrixType M_reg_inv = tf_test_bas2act->GetInverseMatrix();
	Transform::OffsetType O_reg = tf_test_bas2act->GetOffset();

	Transform::MatrixType M_coreg = tf_coreg->GetMatrix();
	Transform::OffsetType O_coreg = tf_coreg->GetOffset();

	Transform::MatrixType M_final = M_reg_inv * M_coreg;
	Transform::OffsetType O_final = M_reg_inv * (O_coreg - O_reg - (o_bas - o_zero)) + (o_act - o_zero);

	tf_coreg_act = Transform::New();
	tf_coreg_act->SetMatrix(M_final);
	tf_coreg_act->SetOffset(O_final);

	Transform::Pointer tf_coreg_act_i = Transform::New();
	tf_coreg_act->GetInverse(tf_coreg_act_i);

	// ---------------------------------------------------------------------------------------
	// TEST TRANSFORMS CONFIGURATION END -----------------------------------------------------
	// ---------------------------------------------------------------------------------------


	// ---------------------------------------------------------------------------------------
	// GOLD TRANSFORMS CONFIGURATION ---------------------------------------------------------
	// ---------------------------------------------------------------------------------------

	// Basal Gold transform (RM->Interictal): tf_gold_rm2bas
	EulerTransform::Pointer tf_gold_rm2bas;
	itk::TransformFileReader::Pointer tf_gold_rm2bas_reader = itk::TransformFileReader::New();
	tf_gold_rm2bas_reader->SetFileName(bas_tfs_file);
	tf_gold_rm2bas_reader->Update();

	TFList tfs_bas_list = tf_gold_rm2bas_reader->GetTransformList();
	itk::TransformFileReader::TransformListType::const_iterator tfs_bas_it = tfs_bas_list->begin();

	std::advance(tfs_bas_it, (case_id - 1));

	if (!strcmp((*tfs_bas_it)->GetNameOfClass(), "Euler3DTransform")) {
		tf_gold_rm2bas = static_cast<EulerTransform *> ((*tfs_bas_it).GetPointer());
	} else {
		std::cout << "Reading basal gold transform failed" << std::endl;
	}

	EulerTransform::Pointer tf_gold_rm2bas_inv = EulerTransform::New();
	tf_gold_rm2bas_inv->SetCenter(tf_gold_rm2bas->GetCenter());
	tf_gold_rm2bas->GetInverse(tf_gold_rm2bas_inv);

	// ---------------------------------------------------------------------------
	// Activated Gold transform (RM->Ictal): tf_act
	EulerTransform::Pointer tf_gold_rm2act;
	itk::TransformFileReader::Pointer tf_gold_rm2act_reader = itk::TransformFileReader::New();
	tf_gold_rm2act_reader->SetFileName(act_tfs_file);
	tf_gold_rm2act_reader->Update();

	TFList tfs_act_list = tf_gold_rm2act_reader->GetTransformList();
	itk::TransformFileReader::TransformListType::const_iterator tfs_act_it = tfs_act_list->begin();

	std::advance(tfs_act_it, (case_id - 1));

	if (!strcmp((*tfs_act_it)->GetNameOfClass(), "Euler3DTransform")) {
		tf_gold_rm2act = static_cast<EulerTransform *> ((*tfs_act_it).GetPointer());
	} else {
		std::cout << "Reading activated gold transform failed" << std::endl;
	}

	// Real space's center Correction
	EulerTransform::Pointer t = EulerTransform::New();
	t->SetIdentity();
	t->SetTranslation(-v_rm2ii);
	tf_gold_rm2bas_inv->Compose(t, false);

	EulerTransform::Pointer tf_gold_rm2act_inv = EulerTransform::New();
	tf_gold_rm2act_inv->SetCenter(tf_gold_rm2act->GetCenter());
	tf_gold_rm2act->GetInverse(tf_gold_rm2act_inv);

	// Real space's center Correction
	EulerTransform::Pointer t2 = EulerTransform::New();
	t2->SetIdentity();
	t2->SetTranslation(-v_rm2i);
	tf_gold_rm2act_inv->Compose(t2, false);

	// ---------------------------------------------------------------------------------------
	// GOLD TRANSFORMS CONFIGURATION END -----------------------------------------------------
	// ---------------------------------------------------------------------------------------


	if (output_all || output_ref) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage(true);
		res_bas->SetReferenceImage(i_ref);
		res_bas->SetInput(i_bas);
		res_bas->SetTransform(tf_coreg);
		res_bas->Update();
		w->SetInput(res_bas->GetOutput());
		w->SetFileName(outputPrefix + "interictal_test.nii.gz");
		w->Update();

		Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage(true);
		res_act->SetReferenceImage(i_ref);
		res_act->SetInput(i_act);
		res_act->SetTransform(tf_coreg_act);
		res_act->Update();
		w->SetInput(res_act->GetOutput());
		w->SetFileName(outputPrefix + "ictal_test.nii.gz");
		w->Update();
	}

	if (output_all || output_gld) {
		Resampler::Pointer res_gold_ii = Resampler::New();
		res_gold_ii->SetInput(i_bas);
		res_gold_ii->SetUseReferenceImage(true);
		res_gold_ii->SetReferenceImage(i_ref);
		res_gold_ii->SetTransform(tf_gold_rm2bas_inv);
		res_gold_ii->Update();
		w->SetInput(res_gold_ii->GetOutput());
		w->SetFileName(outputPrefix + "interictal_gold.nii.gz");
		w->Update();

		Resampler::Pointer res_gold_i = Resampler::New();
		res_gold_i->SetInput(i_act);
		res_gold_i->SetUseReferenceImage(true);
		res_gold_i->SetReferenceImage(i_ref);
		res_gold_i->SetTransform(tf_gold_rm2act_inv);
		res_gold_i->Update();
		w->SetInput(res_gold_i->GetOutput());
		w->SetFileName(outputPrefix + "ictal_gold.nii.gz");
		w->Update();
	}

	if (output_all || output_bas) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage(false);
		res_bas->SetOutputParametersFromImage(i_bas);
		res_bas->SetSize(s_bas);
		res_bas->SetOutputOrigin(o_bas + v_rm2ii);
		res_bas->SetInput(i_bas);
		res_bas->SetTransform(tf_coreg);
		res_bas->Update();
		w->SetInput(res_bas->GetOutput());
		w->SetFileName(outputPrefix + "interictal_bas_test.nii.gz");
		w->Update();

		Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage(false);
		res_act->SetOutputParametersFromImage(i_bas);
		res_act->SetSize(s_bas);
		res_act->SetOutputOrigin(o_bas + v_rm2ii);
		res_act->SetInput(i_act);
		res_act->SetTransform(tf_coreg_act);
		res_act->Update();
		w->SetInput(res_act->GetOutput());
		w->SetFileName(outputPrefix + "ictal_bas_test.nii.gz");
		w->Update();
	}

	if (output_all || output_act) {
		Resampler::Pointer res_bas = Resampler::New();
		res_bas->SetUseReferenceImage(false);
		res_bas->SetOutputParametersFromImage(i_act);
		res_bas->SetSize(s_act);
		res_bas->SetOutputOrigin(o_bas + v_rm2i);
		res_bas->SetInput(i_bas);
		res_bas->SetTransform(tf_coreg);
		res_bas->Update();

		w->SetInput(res_bas->GetOutput());
		w->SetFileName(outputPrefix + "interictal_act_test.nii.gz");
		w->Update();

		Resampler::Pointer res_act = Resampler::New();
		res_act->SetUseReferenceImage(false);
		res_act->SetOutputParametersFromImage(i_act);
		res_act->SetSize(s_act);
		res_act->SetOutputOrigin(o_act + v_rm2i);
		res_act->SetInput(i_act);
		res_act->SetTransform(tf_coreg_act);
		res_act->Update();

		w->SetInput(res_act->GetOutput());
		w->SetFileName(outputPrefix + "ictal_act_test.nii.gz");
		w->Update();
	}

	if (output_all || output_dummy) {
		Resampler::Pointer res_ref_ii = Resampler::New();
		res_ref_ii->SetInput(i_bas);
		res_ref_ii->SetUseReferenceImage(false);
		res_ref_ii->SetOutputParametersFromImage(i_ref);
		res_ref_ii->SetOutputOrigin(bas_ref_orig);
		res_ref_ii->Update();
		w->SetInput(res_ref_ii->GetOutput());
		w->SetFileName(outputPrefix + "bas_ref.nii.gz");
		w->Update();

		Resampler::Pointer res_bas_i = Resampler::New();
		res_bas_i->SetInput(i_act);
		res_bas_i->SetUseReferenceImage(true);
		res_bas_i->SetReferenceImage(i_bas);
		res_bas_i->SetTransform(tf_test_bas2act_inv);
		res_bas_i->Update();
		w->SetInput(res_bas_i->GetOutput());
		w->SetFileName(outputPrefix + "ictal_reg_bas.nii.gz");
		w->Update();

		Resampler::Pointer res_act_ii = Resampler::New();
		res_act_ii->SetInput(i_bas);
		res_act_ii->SetUseReferenceImage(true);
		res_act_ii->SetReferenceImage(i_act);
		res_act_ii->SetTransform(tf_test_bas2act);
		res_act_ii->Update();
		w->SetInput(res_act_ii->GetOutput());
		w->SetFileName(outputPrefix + "interictal_reg_act.nii.gz");
		w->Update();
	}

	//########################################################################################
	// TRE  ----------------------------------------------------------------------------------
	// Berta Martí 2011-03-02

	if (!lesionMask) {
		std::vector<double> errors;
		std::vector<Image::IndexType> cornersIdxs;
		std::vector<Image::IndexType> targetsIdxs;

		//We defined six points
		cornersIdxs.resize(6);
		targetsIdxs.resize(6);

		MomentsCalculator::Pointer getMoment = MomentsCalculator::New();
		getMoment->SetImage(i_brain);
		getMoment->Compute();

		//Get the center of mask of Brain image
		Image::PointType cMass;
		cMass.Fill(0.0);
		cMass += getMoment->GetCenterOfGravity();

		//Generate a binary mask of Brain image
		typedef itk::BinaryThresholdImageFilter<Image, Image> MaskGenerator;
		MaskGenerator::Pointer maskG = MaskGenerator::New();
		maskG->SetInput(i_brain);
		maskG->SetLowerThreshold(1);
		maskG->SetOutsideValue(0);
		maskG->Update();

		//Compute contour image Filter
		typedef itk::BinaryContourImageFilter<Image, Image> Contour;
		Contour::Pointer binaryContourFilter = Contour::New();
		binaryContourFilter->SetInput(maskG->GetOutput());
		binaryContourFilter->Update();

		//Get the TRE points (6, one in each direction)
		Image::Pointer edgeImage = binaryContourFilter->GetOutput();
		edgeImage->SetDirection(i_ref->GetDirection());

		Image::IndexType indexCM;
		edgeImage->TransformPhysicalPointToIndex(cMass, indexCM);

		if (output_msk) {
			mask_out->Allocate();
			mask_out->FillBuffer(0.0);
			mask_out->SetPixel(indexCM, 750);
		}

		for (unsigned int i = 0; i < cornersIdxs.size(); i++) {
			cornersIdxs[i] = indexCM;

			for (unsigned int j = 0; j < Dimension; j++) {
				if (i == j) {
					cornersIdxs[i][j] = edgeImage->GetLargestPossibleRegion().GetSize()[j] - 1;
				}
				if (i == j + Dimension) {
					cornersIdxs[i][j] = 0;
				}
			}

			itk::LineConstIterator<Image> itLine(edgeImage, indexCM, cornersIdxs[i]);
			itk::LineConstIterator<Image> itOut(i_brain, indexCM, cornersIdxs[i]);
			itLine.GoToBegin();
			itOut.GoToBegin();

			bool limit = false;

			while (!itLine.IsAtEnd()) {
				if (itLine.Get() > 0.0) {
					targetsIdxs[i] = itOut.GetIndex();
					limit = true;
				}

				++itLine;
				++itOut;

				//If itLine arrives to the end of the line and no pixel is 255 (limit==false) put it on the corner
				if (itLine.IsAtEnd() && limit == false) {
					targetsIdxs[i] = cornersIdxs[i];
					//To avoid boundary problems, we move two pixels the position on the limit of the image
					for (int j = 0; j < 3; j++) {
						if (targetsIdxs[i][j] == 0) {
							targetsIdxs[i][j] = 2;
						}

						if (targetsIdxs[i][j] == edgeImage->GetLargestPossibleRegion().GetSize()[j]) {
							targetsIdxs[i][j] = edgeImage->GetLargestPossibleRegion().GetSize()[j] - 2;
						}

					}
				}

				//If itLine arrives to the end of the line and the pixel 255 is on the corner put it on the corner + 2 or -2
				if (itLine.IsAtEnd() && limit == true) {
					for (int j = 0; j < 3; j++) {
						if (targetsIdxs[i][j] == 0) {
							targetsIdxs[i][j] = 2;
						}

						if (targetsIdxs[i][j] == edgeImage->GetLargestPossibleRegion().GetSize()[j]) {
							targetsIdxs[i][j] = edgeImage->GetLargestPossibleRegion().GetSize()[j] - 2;
						}
					}
				}
			}
		}

		std::vector<double> regError;
		std::vector<double> coregError;
		std::vector<double> totalError;

		for (unsigned int i = 0; i < targetsIdxs.size(); i++) {
			Image::PointType p_ref;
			i_ref->TransformIndexToPhysicalPoint(targetsIdxs[i], p_ref);

			// Buscar los puntos de lesión en el espacio interictal (a partir de la teorica)
			Image::PointType p_bas_gold = tf_gold_rm2bas_inv->TransformPoint(p_ref);
			// Buscar los puntos de lesión en el espacio Ictal (a partir de la teorica)
			Image::PointType p_act_gold = tf_gold_rm2act_inv->TransformPoint(p_ref);

			// Buscar los puntos de lesión en el espacio Interictal (experimental)
			Image::PointType p_coreg_test_ref = tf_test_bas2rm_inv->TransformPoint(p_bas_gold);

			// Buscar los puntos de la lesión en el espacio Ictal (experimental)
			Image::PointType p_reg_test_act = tf_test_bas2act_inv->TransformPoint(p_bas_gold);

			// Buscar los puntos de la lesión en el espacio referencia (experimental)
			//Image::PointType p_total_test_ref = tf_test_act2rm_inv->TransformPoint( p_act_gold );
			Image::PointType p_total_test_ref = tf_coreg_act_i->TransformPoint(p_act_gold);
			//Image::PointType p_total_test_ref =	tf_test_bas2rm_inv->TransformPoint(	tf_test_bas2act->TransformPoint(p_act_gold));

			if (output_msk)
				mask_out->SetPixel(targetsIdxs[i], 1000);

			regError.push_back(p_act_gold.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_reg_test_act));
			coregError.push_back(p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_coreg_test_ref));
			totalError.push_back(p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_total_test_ref));
		}

		std::cout << std::setw(5) << std::left << case_id;
		std::cout << std::setw(10) << std::right << mean(regError);
		std::cout << std::setw(10) << std::right << median(regError);
		std::cout << std::setw(10) << std::right << max(regError) << std::setw(3) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << mean(coregError);
		std::cout << std::setw(10) << std::right << median(coregError);
		std::cout << std::setw(10) << std::right << max(coregError) << std::setw(3) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << mean(totalError);
		std::cout << std::setw(10) << std::right << median(totalError);
		std::cout << std::setw(10) << std::right << max(totalError) << std::setw(3) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << totalError.size() << std::endl;

		if (output_msk) {
			typedef itk::ConstantBoundaryCondition<Image> BoundaryConditionType;
			typedef itk::NeighborhoodIterator<Image, BoundaryConditionType> NeighborhoodIteratorType;
			//Fill the six points
			NeighborhoodIteratorType::RadiusType radius;
			radius.Fill(2);
			NeighborhoodIteratorType itNB(radius, mask_out, mask_out->GetLargestPossibleRegion());
			std::vector<NeighborhoodIteratorType::OffsetType> offsetIt;
			offsetIt.resize(125);

			BoundaryConditionType boundaryCondition;
			boundaryCondition.SetConstant(0);
			itNB.SetBoundaryCondition(boundaryCondition);

			int iniOffsetIt = 0;

			//Fill offset vector in order to have all neighborhood offsets (27)
			for (unsigned int i = -2; i < Dimension; i++) {
				for (unsigned int j = -2; j < Dimension; j++) {
					for (unsigned int k = -2; k < Dimension; k++) {
						offsetIt[iniOffsetIt][0] = i;
						offsetIt[iniOffsetIt][1] = j;
						offsetIt[iniOffsetIt][2] = k;
						iniOffsetIt++;
					}
				}
			}

			//If mean value of neighborhood values is four times larger than centered value, put this one zero.
			for (itNB.GoToBegin(); !itNB.IsAtEnd(); ++itNB) {
				vnl_vector<float> vectorValues;
				vectorValues.set_size(offsetIt.size());

				if (itNB.GetCenterPixel() == 1000) {
					for (unsigned int i = 0; i < offsetIt.size(); i++) {
						itNB.SetPixel(offsetIt[i], 500);
					}
				}

				if (itNB.GetCenterPixel() == 750) {
					for (unsigned int i = 0; i < offsetIt.size(); i++) {
						itNB.SetPixel(offsetIt[i], 600);
					}
				}
			}

			mask_out->Update();
			w->SetInput(mask_out);
			w->SetFileName(outputPrefix + "MaskTRE_0" + s_id.str() + ".nii.gz");
			w->Update();
		}

	}
	//---------------------------------------------------------------------------------------

	if (lesionMask) {
		Image::PointType cm_ref; cm_ref.Fill(0.0);

		i_mask->SetOrigin( cm_ref );

		MomentsCalculator::Pointer calc_les = MomentsCalculator::New();
		calc_les->SetImage(i_mask);
		calc_les->Compute();
		cm_ref += calc_les->GetCenterOfGravity();

		double norm = center_ref.EuclideanDistanceTo<Image::PointType::CoordRepType>(cm_ref);

		std::vector<double> regError;
		std::vector<double> coregError;
		std::vector<double> totalError;

		Iterator it(i_mask, i_mask->GetLargestPossibleRegion());

		if (output_msk) {
			mask_out->Allocate();
			mask_out->FillBuffer(0.0);
		}

		Image::PointType cm_test; cm_test.Fill(0.0);

		for (it = it.Begin(); !it.IsAtEnd(); ++it) {
			if (it.Get() > 0.0) {
				Image::PointType p_ref;
				i_ref->TransformIndexToPhysicalPoint(it.GetIndex(), p_ref);

				cm_test[0]+= p_ref[0];
				cm_test[1]+= p_ref[1];
				cm_test[2]+= p_ref[2];

				// Buscar los puntos de lesión en el espacio interictal (a partir de la teorica)
				Image::PointType p_bas_gold = tf_gold_rm2bas_inv->TransformPoint(p_ref);
				// Buscar los puntos de lesión en el espacio Ictal (a partir de la teorica)
				Image::PointType p_act_gold = tf_gold_rm2act_inv->TransformPoint(p_ref);

				// Buscar los puntos de lesión en el espacio Interictal (experimental)
				Image::PointType p_coreg_test_ref = tf_test_bas2rm_inv->TransformPoint(p_bas_gold);

				// Buscar los puntos de la lesión en el espacio Ictal (experimental)
				Image::PointType p_reg_test_act = tf_test_bas2act_inv->TransformPoint(p_bas_gold);

				// Buscar los puntos de la lesión en el espacio referencia (experimental)
				//Image::PointType p_total_test_ref = tf_test_act2rm_inv->TransformPoint( p_act_gold );
				Image::PointType p_total_test_ref = tf_test_bas2rm_inv->TransformPoint(tf_test_bas2act->TransformPoint(p_act_gold));

				regError.push_back(p_act_gold.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_reg_test_act));
				coregError.push_back(p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_coreg_test_ref));
				totalError.push_back(p_ref.EuclideanDistanceTo<Image::PointType::CoordRepType> (p_total_test_ref));

				if (output_msk)
					mask_out->SetPixel(it.GetIndex(), 1);
			}
		}

		cm_test[0] = cm_test[0] / regError.size();	
		cm_test[1] = cm_test[1] / regError.size();
		cm_test[2] = cm_test[2] / regError.size();

		std::cout << std::setw(5) << std::left << case_id;
		std::cout << std::setw(10) << std::right << mean(regError);
		std::cout << std::setw(10) << std::right << median(regError);
		std::cout << std::setw(10) << std::right << max(regError) << std::setw(3) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << mean(coregError);
		std::cout << std::setw(10) << std::right << median(coregError);
		std::cout << std::setw(10) << std::right << max(coregError) << std::setw(3) << "|" << std::setw(3) << " ";
		std::cout << std::setw(10) << std::right << mean(totalError);
		std::cout << std::setw(10) << std::right << median(totalError);
		std::cout << std::setw(10) << std::right << max(totalError) << std::setw(2) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << totalError.size() << std::setw(2) << "|" << std::setw(2) << " ";
		std::cout << std::setw(10) << std::right << cm_ref << ", norm=" << norm << std::endl;

		if (output_msk) {
			mask_out->Update();
			w->SetInput(mask_out);
			w->SetFileName(outputPrefix + "MaskLRE0" + s_id.str() + ".nii.gz");
			w->Update();
		}
	}

	return EXIT_SUCCESS;
}
