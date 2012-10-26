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

/*
 * cumulativeTransform.cxx
 *
 *  Created on: 18/05/2011
 *      Author: oesteban
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

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#define DEG M_PI/180.0

namespace bfs = boost::filesystem;
namespace bpo = boost::program_options;

int main( int argc, char *argv[] )
{
    typedef  unsigned char  PixelRange;
    const int Dimension = 3;

    typedef itk::Image< PixelRange, Dimension >          Image;
    typedef itk::ImageRegionIterator< Image >            Iterator;
    typedef itk::Rigid3DTransform< double >              RigidTransform;
    typedef itk::Euler3DTransform< double >              EulerTransform;
    typedef itk::VersorRigid3DTransform< double >        VersorsTransform;
    typedef itk::ImageFileReader< Image  >               Reader;
    typedef itk::ImageFileWriter< Image >                Writer;
    typedef itk::ResampleImageFilter< Image, Image>      Resampler;

    unsigned int case_id;
    std::string inputFileName;
    std::string outputFileName;
    std::string transformsFileName;
    std::string referenceImage;
    std::string outputPrefix;
    bool doInverse, doFix, doIso = false;
    unsigned int nX, nY, nZ=0;

	bpo::options_description desc("Usage");
	desc.add_options()
			("help", "show help message")
			("input-image", bpo::value< std::string >(&inputFileName)->required(), "input image file")
			("output-image", bpo::value< std::string >(&outputFileName)->required(), "output image file")
			("transforms-file,F", bpo::value< std::string >(&transformsFileName)->required(), "file containing the transforms list")
			("invert,I", bpo::bool_switch( &doInverse )->default_value(false), "invert all transforms")
			("case-id", bpo::value< unsigned int >(&case_id)->default_value(0u), "case list identifier")
			("reference-image,R", bpo::value< std::string >(&referenceImage), "reference image file")
			("fix-volume,V", bpo::bool_switch( &doFix )->default_value(false), "enlarge & resize target volume to be isotropic and isometric")
			("isotropic-output,S", bpo::bool_switch( &doIso )->default_value(false), "make output isotropic. It automatically resizes for the minimun required size")
			("x-val,x", bpo::value< unsigned int >(&nX)->default_value(0u), "specify x size")
			("y-val,y", bpo::value< unsigned int >(&nY)->default_value(0u), "specify y size")
			("z-val,z", bpo::value< unsigned int >(&nZ)->default_value(0u), "specify z size")
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

    std::advance( last_case, case_id );

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

    //Reference
    Image::SizeType ref_size = i_ref->GetLargestPossibleRegion().GetSize();
    Image::SpacingType ref_sp = i_ref->GetSpacing();
    Image::PointType ref_center;   
    double ref_center_idx[3] = { (ref_size[0]-1)*0.5, (ref_size[1]-1)*0.5, (ref_size[2]-1)*0.5 };
    i_ref->TransformContinuousIndexToPhysicalPoint<double> (ref_center_idx, ref_center );

    //initial size can be customized by each dimension
    Image::SizeType out_size; 
    Image::SpacingType out_sp;
    out_size=ref_size;
    if (nX>0){ out_size[0]=nX; doFix=false;}
    if (nY>0){ out_size[1]=nY; doFix=false;}	
    if (nZ>0){ out_size[2]=nZ; doFix=false;}    

    Resampler::Pointer res_gold = Resampler::New();
    //res_gold->SetUseReferenceImage( true );
    if ( doIso ) {		

		unsigned int maxSize = ref_size[0];
		double minSp = ref_sp[0];		
		int ratioSp[3] = {1,1,1};
		int outDim[3] = {0,0,0}; outDim[0]=ref_size[0];outDim[1]=ref_size[1];outDim[2]=ref_size[2];
		for (unsigned int i = 1; i<Dimension; i++) {
				if ( ref_size[i] > maxSize ) maxSize = ref_size[i];
				if ( ref_sp[i] < minSp  ) minSp = ref_sp[i];											
		}

		for (unsigned int i = 1; i<Dimension; i++) {
				if ( ref_sp[i] > minSp  ) {
					ratioSp[i]= (int)((ref_sp[i]/minSp)+0.5);
					std::cout<<"Spacing ratio for dimension "<<i<<" is "<<ratioSp[i]<<std::endl; 	
					outDim[i]=ref_size[i]*ratioSp[i];
				}											
		}

		//We set the minimum as the size required to be isotropic
		if(out_size[0]<outDim[0]) out_size[0]=outDim[0];
		if(out_size[1]<outDim[1]) out_size[1]=outDim[1];
		if(out_size[2]<outDim[2]) out_size[2]=outDim[2];
		out_sp.Fill( minSp );
	    
    }else{
	    	
	    if ( !doFix ) {
	        //no modifications in spacing
	    	out_sp=ref_sp;//res_gold->SetReferenceImage( i_ref );
	    } else {
	    
	 		unsigned int maxSize = ref_size[0];
			double minSp = ref_sp[0];

			for (unsigned int i = 1; i<Dimension; i++) {
				if ( ref_size[i] > maxSize ) maxSize = ref_size[i];
				if ( ref_sp[i] < minSp  ) minSp = ref_sp[i];
			}

			Image::SizeType out_size; out_size.Fill( maxSize );
			Image::SpacingType out_sp; out_sp.Fill( minSp );
			
	    }
    }
    
    Image::PointType out_orig;
    for (unsigned int i = 0; i<Dimension; i++)
        out_orig[i] = ref_center[i] - (out_size[i] * 0.5 * out_sp[i]);
    std::cerr<<"Output image size "<<out_size[0]<<" "<<out_size[1]<<" "<<out_size[2]<<std::endl;
    out_orig = i_ref->GetDirection() * out_orig;
    res_gold->SetOutputDirection( i_ref->GetDirection() );
    res_gold->SetOutputSpacing( out_sp );
    res_gold->SetOutputOrigin( out_orig );
    res_gold->SetSize( out_size );


    res_gold->SetInput( i_in );
    res_gold->SetTransform( tf_gold );
    res_gold->Update();

    w->SetInput( res_gold->GetOutput() );
    w->SetFileName( outputFileName );
    w->Update();

    return EXIT_SUCCESS;
}


