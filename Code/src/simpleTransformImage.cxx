/* --------------------------------------------------------------------------------------
 * File:    simpleTransformImage.cxx
 * Date:    Nov 11, 2011
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2011, Oscar Esteban - oesteban@die.upm.es
 with Biomedical Image Technology, UPM (BIT-UPM)
 and Signal Processing Lab 5, EPFL (LTS5-EPFL)
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

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>
#include <itkRigid3DTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkContinuousIndex.h>
#include <itkImageRegionIterator.h>
#include <itkTransformFileWriter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0

namespace bfs = boost::filesystem;
namespace bpo = boost::program_options;

int main( int argc, char *argv[] )
{
    typedef short  PixelRange;
    const int Dimension = 3;

    typedef itk::Image< PixelRange, Dimension >          Image;
    typedef itk::ImageRegionIterator< Image >            Iterator;
    typedef itk::Rigid3DTransform< double >              RigidTransform;
    typedef itk::Euler3DTransform< double >              EulerTransform;
    typedef itk::VersorRigid3DTransform< double >        VersorsTransform;
    typedef itk::ImageFileReader< Image  >               Reader;
    typedef itk::ImageFileWriter< Image >                Writer;
    typedef itk::ResampleImageFilter< Image, Image>      Resampler;

    unsigned int case_id = 0;
    std::string inputFileName;
    std::string outputFileName;
    std::string transformsFileName;
    std::string referenceImage;
    std::string referenceMovingImage;
    std::string outputPrefix;
    bool doInverse = false;
    bool doFix = false;
    bool useNN = false;
    bool useBSplines = false;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("input-image", bpo::value< std::string >(&inputFileName)->required(), "input image file")
			("output-image", bpo::value< std::string >(&outputFileName)->required(), "output image file")
			("transforms-file,F", bpo::value< std::string >(&transformsFileName)->required(), "file containing the transforms list")
			("invert,I", bpo::bool_switch( &doInverse )->default_value(false), "invert all transforms")
			("use-NN", bpo::bool_switch( &useNN )->default_value(false), "use NN interpolation")
			("use-BSplines", bpo::bool_switch( &useBSplines )->default_value(false), "use b-splines interpolation")
			("reference-image,R", bpo::value< std::string >(&referenceImage), "reference image file")
			("reference-moving,M", bpo::value< std::string >(&referenceMovingImage), "reference moving image file")
			("output-prefix,o", bpo::value< std::string >(&outputPrefix), "prefix for output files");

	bpo::positional_options_description pdesc;
	pdesc.add("input-image", 1);
	pdesc.add("output-image", 2);

	bpo::variables_map vmap;
	bpo::store(bpo::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vmap);

	if (vmap.count("help")) {
	    std::cout << desc << std::endl;
	    return 1;
	}

	try {
		bpo::notify(vmap);
	} catch ( boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::program_options::required_option> > &err){
		std::cout << "Error: " << err.what() << std::endl;
		std::cout << desc << std::endl;
		return EXIT_FAILURE;
	}


    if (!bfs::exists( transformsFileName ) )
    {
      std::cout << "Error: Transforms file not found (supplied path=" << transformsFileName << ")." << std::endl;
      return EXIT_FAILURE;
    }


    // ----------------------------------------------------------------
    // LECTURA DE LAS IMÁGENES
    // ----------------------------------------------------------------
    Reader::Pointer r  = Reader::New();
    Writer::Pointer w  = Writer::New();


    r->SetFileName(  inputFileName );
    try {
        r->Update();
    }
    catch ( itk::ExceptionObject & err )
    {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return EXIT_FAILURE;
    }

    Image::Pointer i_in = r->GetOutput();


    // TRANSFORMS CONFIGURATION ---------------------------------------------------------
    RigidTransform::Pointer tf_gold = RigidTransform::New();
    tf_gold->SetIdentity();

    itk::TransformFileReader::Pointer tf_gold_reader = itk::TransformFileReader::New();
    tf_gold_reader->SetFileName( transformsFileName );
    tf_gold_reader->Update();

    typedef itk::TransformFileReader::TransformListType * TFList;
    TFList tfs_gold_list = tf_gold_reader->GetTransformList();
    itk::TransformFileReader::TransformListType::const_iterator it_gold   = tfs_gold_list->begin();
    itk::TransformFileReader::TransformListType::const_iterator last_case = tfs_gold_list->begin();

    std::advance( last_case, case_id +1 );

    while( it_gold != last_case && it_gold != tfs_gold_list->end() ) {
    	RigidTransform::Pointer tf_step;
    	RigidTransform::Pointer tf_step_inv;

		if (!strcmp( (*it_gold)->GetNameOfClass(), "VersorRigid3DTransform" ))
		{
			tf_step = static_cast< VersorsTransform * > ( (*it_gold).GetPointer() ) ;
		}

		if (doInverse) {
			tf_step_inv = RigidTransform::New();
			tf_step->GetInverse( tf_step_inv );
			tf_step = tf_step_inv;
		}

		tf_gold->Compose( tf_step );

		++it_gold;
    }

    Image::Pointer i_ref = i_in;

    if( bfs::exists( referenceImage ) ) {
    	Reader::Pointer r_ref = Reader::New();
    	r_ref->SetFileName( referenceImage );

    	try {
    		r_ref->Update();
    		i_ref = r_ref->GetOutput();
    	} catch ( itk::ExceptionObject & err ) {
    		std::cout << "Warning: Exception caught (" << err << ")" << std::endl;
    	}
    }

    if( bfs::exists( referenceMovingImage ) ) {
        Image::Pointer i_ref_mov;
    	Reader::Pointer r_ref = Reader::New();
    	r_ref->SetFileName( referenceMovingImage );

    	try {
    		r_ref->Update();
    		i_ref_mov = r_ref->GetOutput();
    	} catch ( itk::ExceptionObject & err ) {
    		std::cout << "Warning: Exception caught (" << err << ")" << std::endl;
    	}
    	double centerIdx[3];

    	Image::SizeType sizeRef = i_ref->GetLargestPossibleRegion().GetSize();
    	Image::PointType centerRef;
    	centerIdx[0] = (sizeRef[0]-1)*0.5; centerIdx[1] = (sizeRef[1]-1)*0.5; centerIdx[2] = (sizeRef[2]-1)*0.5;
    	i_ref->TransformContinuousIndexToPhysicalPoint<Image::PointType::CoordRepType>(centerIdx, centerRef );

    	std::cout << "* Ref image: "    << std::endl;
    	std::cout << "\tSize= "         << sizeRef << std::endl;
    	std::cout << "\tSpacing= "      << i_ref->GetSpacing() << std::endl;
    	std::cout << "\tOrigin= "       << i_ref->GetOrigin() << std::endl;
    	std::cout << "\tCenter= "       << centerRef << std::endl;
    	std::cout << "\tDirection= "    << std::endl << i_ref->GetDirection() << std::endl << std::endl;


        Image::SizeType sizeMov = i_in->GetLargestPossibleRegion().GetSize();
        Image::SpacingType spMov = i_in->GetSpacing();
        Image::PointType oMov = i_in->GetOrigin();

        Image::PointType centerMov;
        centerIdx[0] = (sizeMov[0]-1)*0.5; centerIdx[1] = (sizeMov[1]-1)*0.5; centerIdx[2] = (sizeMov[2]-1)*0.5;
        i_in->TransformContinuousIndexToPhysicalPoint<Image::PointType::CoordRepType>(centerIdx, centerMov );

        std::cout << "* Moving image: " << std::endl;
        std::cout << "\tSize= "         << sizeMov << std::endl;
        std::cout << "\tSpacing= "      << spMov << std::endl;
        std::cout << "\tOrigin= "       << oMov << std::endl;
        std::cout << "\tCenter= "       << centerMov << std::endl;
        std::cout << "\tDirection= "    << std::endl << i_in->GetDirection() << std::endl << std::endl;


        Image::SizeType sizeRefMov = i_ref_mov->GetLargestPossibleRegion().GetSize();
        Image::SpacingType spRefMov = i_ref_mov->GetSpacing();
        Image::PointType oRefMov = i_ref_mov->GetOrigin();

        Image::PointType centerRefMov;
        centerIdx[0] = (sizeRefMov[0]-1)*0.5; centerIdx[1] = (sizeRefMov[1]-1)*0.5; centerIdx[2] = (sizeRefMov[2]-1)*0.5;
        i_ref_mov->TransformContinuousIndexToPhysicalPoint<Image::PointType::CoordRepType>(centerIdx, centerRefMov );

        std::cout << "* RefMov image: " << std::endl;
        std::cout << "\tSize= "         << sizeRefMov << std::endl;
        std::cout << "\tSpacing= "      << spRefMov << std::endl;
        std::cout << "\tOrigin= "       << oRefMov << std::endl;
        std::cout << "\tCenter= "       << centerRefMov << std::endl;
        std::cout << "\tDirection= "    << std::endl << i_ref_mov->GetDirection() << std::endl << std::endl;


        //VersorsTransform::InputVectorType t = centerMov - centerRefMov;
        //VersorsTransform::InputVectorType t2 = centerMov - centerRefMov;
        VersorsTransform::InputVectorType t2 = centerRefMov - centerMov;
        //Image::PointType newOriginMov = oMov - t;

        VersorsTransform::InputVectorType t_tfed = tf_gold->TransformVector( t2 );
        VersorsTransform::InputVectorType t3 = centerRef - tf_gold->TransformPoint( centerMov );

        VersorsTransform::Pointer tf = VersorsTransform::New();
        tf->SetIdentity();
        tf->SetOffset( t2 );

        tf_gold->Compose( tf, true );


        VersorsTransform::Pointer tf_post = VersorsTransform::New();
        tf_post->SetIdentity();
        tf_post->SetOffset( - t_tfed );
        //tf_post->SetOffset( t3 );

        VersorsTransform::Pointer tf_inv = VersorsTransform::New();
        tf_post->GetInverse( tf_inv );

        //tf_gold->Compose( tf_inv );
        tf_gold->Compose( tf_post );

    }


    Resampler::Pointer res_gold = Resampler::New();
    res_gold->SetUseReferenceImage( true );

    if ( useNN ) {
    	res_gold->SetInterpolator( itk::NearestNeighborInterpolateImageFunction< Image >::New() );
    } else if ( useBSplines ) {
    	res_gold->SetInterpolator( itk::BSplineInterpolateImageFunction< Image >::New() );
	}

    if ( !doFix ) {
    	res_gold->SetReferenceImage( i_ref );
    } else {
    	Image::SizeType ref_size = i_ref->GetLargestPossibleRegion().GetSize();
    	Image::SpacingType ref_sp = i_ref->GetSpacing();
    	Image::PointType ref_center;

    	double ref_center_idx[3] = { (ref_size[0]-1)*0.5, (ref_size[1]-1)*0.5, (ref_size[2]-1)*0.5 };

		i_ref->TransformContinuousIndexToPhysicalPoint<double> (ref_center_idx, ref_center );

		unsigned int maxSize = ref_size[0];
		double minSp = ref_sp[0];

		for (unsigned int i = 1; i<Dimension; i++) {
			if ( ref_size[i] > maxSize ) maxSize = ref_size[i];
			if ( ref_sp[i] < minSp  ) minSp = ref_sp[i];
		}

		Image::SizeType out_size; out_size.Fill( maxSize );
		Image::SpacingType out_sp; out_sp.Fill( minSp );

		Image::PointType out_orig;

		for (unsigned int i = 0; i<Dimension; i++)
			out_orig[i] = ref_center[i] - (out_size[i] * 0.5 * out_sp[i]);

		out_orig = i_ref->GetDirection() * out_orig;

		res_gold->SetOutputDirection( i_ref->GetDirection() );
		res_gold->SetOutputSpacing( out_sp );
		res_gold->SetOutputOrigin( out_orig );
		res_gold->SetSize( out_size );
    }
    res_gold->SetInput( i_in );
    res_gold->SetTransform( tf_gold );
    res_gold->Update();

    w->SetInput( res_gold->GetOutput() );
    w->SetFileName( outputFileName );
    w->Update();

    return EXIT_SUCCESS;
}



