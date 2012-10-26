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

// -------------------------------------------------------------------------------
// transformITKImage.cxx
// v1.0 - 2011-01-27
//
// This code performs a geometric transform on the input image to an output image
// using the transform number case_id on the file transformsItkFile.
// This file must have the ITKTransform list format.
// -------------------------------------------------------------------------------

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkVersorRigid3DTransform.h>
#include <itkGIBUB3DTransform.h>
#include <itkOrientedImage.h>
#include <itkOrientImageFilter.h>
#include <itkMaskImageFilter.h>

#include <itkContinuousIndex.h>
#include <itkImageRegionIterator.h>
#include <itkTransformFileWriter.h>
#include <itkImageMaskSpatialObject.h>

#include "itkSISCOMSubtractionFilter.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0


namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;


int main( int argc, char *argv[] )
{
    typedef float SpectPixel;
    typedef char           DiffPixel;
    typedef short          MaskPixel;
    const int Dimension = 3;

    typedef itk::Image< SpectPixel, Dimension >                SpectImg;
    typedef itk::Image< DiffPixel, Dimension >                 DiffImg;
    typedef itk::Image< MaskPixel, Dimension >                 Mask;
    typedef itk::ImageRegionIterator< SpectImg >               Iterator;
    typedef itk::VersorRigid3DTransform< double >              Transform;
    typedef itk::ImageFileReader< SpectImg  >                  SpectReader;
    typedef itk::ImageFileReader< DiffImg  >                   DiffReader;
    typedef itk::ImageFileReader< Mask  >                      MaskReader;
    typedef itk::ImageFileWriter< SpectImg >                   SpectWriter;
    typedef itk::ImageFileWriter< DiffImg >                    DiffWriter;
    typedef itk::ResampleImageFilter< SpectImg, SpectImg>      Resampler;
    typedef itk::ResampleImageFilter< DiffImg, DiffImg>        DiffResampler;
    typedef itk::OrientImageFilter< SpectImg, SpectImg >       Orienter;

    typedef itk::MaskImageFilter< DiffImg, Mask, DiffImg >     DiffMasker;


    std::string interictalFileName, interictalMaskName;
    std::string ictalFileName, ictalMaskName;
    std::string mskFileName;
    std::string referenceFileName;
    std::string transformsFileName;
    std::string outputPrefix;
    bool output_images = false;
    bool output_ref = true;
    bool output_std = false;
    bool output_bas = false;
    bool output_act = false;
    Mask::Pointer i_msk;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("interictal-image", bpo::value< std::string >(&interictalFileName)->required(), "interictal image file (UChar)")
			("ictal-image", bpo::value< std::string >(&ictalFileName)->required(), "ictal image file (UChar)")
			("interictal-mask", bpo::value< std::string >(&interictalMaskName), "interictal mask (UChar)")
			("ictal-mask", bpo::value< std::string >(&ictalMaskName), "ictal mask (UChar)")
			("reference-image,R", bpo::value< std::string >(&referenceFileName)->required(), "reference image file")
			("reference-mask,M", bpo::value< std::string>(&mskFileName), "mask (in reference space) file")
			("transforms-file,F"
					"", bpo::value< std::string >(&transformsFileName)->required(), "itk::VersorRigid3DTranform transforms file in ITK format. First transform=coregistration. Second transform=registration.")
			("output-prefix,o", bpo::value< std::string >(&outputPrefix), "prefix for output files")
			("output-all,a", bpo::bool_switch(&output_images)->default_value(false), "output all images")
			("output-ref,r", bpo::bool_switch(&output_ref)->default_value(true), "output images in reference space")
			("output-bas,b", bpo::bool_switch(&output_bas)->default_value(false), "output images in basal space")
			("output-act,i", bpo::bool_switch(&output_act)->default_value(false), "output images in activated space")
			("output-std,s", bpo::bool_switch(&output_std)->default_value(false), "output images in standard space [128 128 54] [3.32 3.32 3.32]");

	bpo::positional_options_description pdesc;
	pdesc.add("interictal-image", 1);
	pdesc.add("ictal-image", 2);

	bpo::variables_map vmap;
	bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vmap);

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

	if (bfs::exists(mskFileName)) {
		MaskReader::Pointer r_msk = MaskReader::New();
		r_msk->SetFileName(mskFileName);
		i_msk = r_msk->GetOutput();

		try {
			i_msk->Update();
		} catch (...) {

		}
	}

	Mask::Pointer i_msk_I;
	if (bfs::exists(ictalMaskName)) {
		MaskReader::Pointer r_msk = MaskReader::New();
		r_msk->SetFileName(ictalMaskName);
		i_msk_I = r_msk->GetOutput();

		try {
			i_msk_I->Update();
		} catch (...) {

		}
	}

	Mask::Pointer i_msk_II;
	if (bfs::exists(interictalMaskName)) {
		MaskReader::Pointer r_msk = MaskReader::New();
		r_msk->SetFileName(ictalMaskName);
		i_msk_II = r_msk->GetOutput();

		try {
			i_msk_II->Update();
		} catch (...) {

		}
	}

	if (!bfs::exists(transformsFileName)) {
		std::cout << "Error: Transforms file not found (supplied path=" << transformsFileName << ")." << std::endl;
		return EXIT_FAILURE;
	}

	// ----------------------------------------------------------------
	// LECTURA DE LAS IMÁGENES
	// ----------------------------------------------------------------
	SpectReader::Pointer r_bas = SpectReader::New();
	SpectReader::Pointer r_act = SpectReader::New();
	SpectReader::Pointer r_ref = SpectReader::New();
	SpectWriter::Pointer w = SpectWriter::New();
	DiffWriter::Pointer w_diff = DiffWriter::New();

	r_bas->SetFileName(interictalFileName);
	r_act->SetFileName(ictalFileName);
	r_ref->SetFileName(referenceFileName);
	try {
		r_bas->Update();
		r_act->Update();
		r_ref->Update();
	} catch (itk::ExceptionObject & err) {
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	SpectImg::Pointer i_bas = r_bas->GetOutput();
	SpectImg::Pointer i_act = r_act->GetOutput();
	SpectImg::Pointer i_ref = r_ref->GetOutput();

	SpectImg::PointType o_bas = i_bas->GetOrigin();
	SpectImg::PointType o_act = i_act->GetOrigin();
	SpectImg::PointType o_ref = i_ref->GetOrigin();

	SpectImg::SizeType s_bas = i_bas->GetLargestPossibleRegion().GetSize();
	SpectImg::PointType center_bas;
	double center_index_bas[3];
	center_index_bas[0] = (s_bas[0] - 1) * 0.5;
	center_index_bas[1] = (s_bas[1] - 1) * 0.5;
	center_index_bas[2] = (s_bas[2] - 1) * 0.5;
	i_bas->TransformContinuousIndexToPhysicalPoint<double> (center_index_bas, center_bas);

	SpectImg::SizeType s_ref = i_ref->GetLargestPossibleRegion().GetSize();
	SpectImg::PointType center_ref;
	double center_index_ref[3];
	center_index_ref[0] = (s_ref[0] - 1) * 0.5;
	center_index_ref[1] = (s_ref[1] - 1) * 0.5;
	center_index_ref[2] = (s_ref[2] - 1) * 0.5;
	i_ref->TransformContinuousIndexToPhysicalPoint<double> (center_index_ref, center_ref);

	SpectImg::SizeType s_act = i_act->GetLargestPossibleRegion().GetSize();
	SpectImg::PointType center_act;
	double center_index_act[3];
	center_index_act[0] = (s_act[0] - 1) * 0.5;
	center_index_act[1] = (s_act[1] - 1) * 0.5;
	center_index_act[2] = (s_act[2] - 1) * 0.5;
	i_act->TransformContinuousIndexToPhysicalPoint<double> (center_index_act, center_act);

	SpectImg::PointType o_zero;
	o_zero.Fill(0.0);

	Transform::InputVectorType ov_ii2i = o_bas - o_act;
	Transform::InputVectorType ov_rm2i = o_ref - o_act;
	Transform::InputVectorType ov_rm2ii = o_ref - o_bas;

	Transform::InputVectorType v_rm2ii = center_ref - center_bas;
	Transform::InputVectorType v_rm2i = center_ref - center_act;
	Transform::InputVectorType v_ii2i = center_bas - center_act;

	// TRANSFORMS CONFIGURATION ---------------------------------------------------------
	Transform::Pointer tf_coreg; // Transform  ii 2 RM (co-registration)
	Transform::Pointer tf_reg; // Transform: ii 2 i  (registration)
	Transform::Pointer tf_coreg_i; // Transform: i  2 RM (co-registration combined with registration)

	itk::TransformFileReader::Pointer tf_gold_reader = itk::TransformFileReader::New();
	tf_gold_reader->SetFileName(transformsFileName);
	tf_gold_reader->Update();

	typedef itk::TransformFileReader::TransformListType * TFList;
	TFList tfs_gold_list = tf_gold_reader->GetTransformList();
	itk::TransformFileReader::TransformListType::const_iterator it_gold = tfs_gold_list->begin();

	Transform::Pointer tf_reg_inv = Transform::New();
	tf_reg_inv->SetIdentity();

	try {
		if (!strcmp((*it_gold)->GetNameOfClass(), "VersorRigid3DTransform")) {
			tf_coreg = static_cast<Transform *> ((*it_gold).GetPointer());
		} else {
			std::cout << "Error: No co-registration transform found" << std::endl;
		}

		it_gold++;
	} catch (itk::ExceptionObject & err) {

		std::cout << "Error reading transforms" << std::endl << err << std::endl;
	}

	SpectImg::PointType new_origin;
	new_origin.Fill(0.0);
	new_origin[1] = (i_act->GetLargestPossibleRegion().GetSize()[1] - 1) * i_act->GetSpacing()[1];
	Transform::InputVectorType v = i_act->GetOrigin() - new_origin;
	tf_reg_inv->SetCenter(tf_reg_inv->GetCenter() + v);


	typedef itk::ImageMaskSpatialObject<3> MaskSpatialObjectType;
	typedef itk::CastImageFilter<Mask, MaskSpatialObjectType::ImageType> MaskObjectCaster;

	typedef itk::SISCOMSubtractionFilter<SpectImg, DiffImg> SISCOMFilter;
	SISCOMFilter::Pointer siscom = SISCOMFilter::New();
	siscom->SetIctalImage(i_act);
	siscom->SetInterictalImage(i_bas);
	siscom->SetRegTransform(tf_reg_inv);

	if (i_msk_I.IsNotNull() ) {
		MaskObjectCaster::Pointer cast = MaskObjectCaster::New();
		cast->SetInput(i_msk_I);
		cast->Update();
		MaskSpatialObjectType::Pointer ictal_mask = MaskSpatialObjectType::New();
		ictal_mask->SetImage(cast->GetOutput());
		siscom->SetIctalImageMask(ictal_mask);
	}

	if (i_msk_II.IsNotNull() ) {
		MaskObjectCaster::Pointer cast = MaskObjectCaster::New();
		cast->SetInput(i_msk_II);
		cast->Update();
		MaskSpatialObjectType::Pointer interictal_mask = MaskSpatialObjectType::New();
		interictal_mask->SetImage(cast->GetOutput());
		siscom->SetInterictalImageMask(interictal_mask);
	}

	siscom->SetInterictalTransform( tf_coreg );
	siscom->Update();

	if (output_images || output_act ) {
		w->SetInput(siscom->GetInterictalNormalizedOutput());
		w->SetFileName(outputPrefix + "inter_norm_act.nii.gz");
		w->Update();

		w_diff->SetInput(siscom->GetOutput());
		w_diff->SetFileName(outputPrefix + "siscom_act.nii.gz");
		w_diff->Update();
	}

	if (output_images || output_ref) {
		siscom->ResampleToReferenceSpace(i_ref);

		SISCOMFilter::InputImageConstPointer i_bas_out = siscom->GetInterictalRefSpace();
		SISCOMFilter::InputImageConstPointer i_act_out = siscom->GetIctalRefSpace();
		SISCOMFilter::OutputImageConstPointer i_dif_out  = siscom->GetSubtractionRefSpace();

		w->SetInput(i_bas_out);
		w->SetFileName(outputPrefix + "inter_norm_ref.nii.gz");
		w->Update();

		w->SetInput(i_act_out);
		w->SetFileName(outputPrefix + "ictal_ref.nii.gz");
		w->Update();

		if (i_msk.IsNotNull()) {
			DiffMasker::Pointer s_masker = DiffMasker::New();
			s_masker->SetInput1(i_dif_out);
			s_masker->SetInput2(i_msk);
			s_masker->Update();
			i_dif_out = s_masker->GetOutput();
		}
		w_diff->SetInput(i_dif_out);
		w_diff->SetFileName(outputPrefix + "siscom_ref.nii.gz");
		w_diff->Update();

	}

	return EXIT_SUCCESS;
}
