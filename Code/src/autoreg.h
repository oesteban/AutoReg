/* --------------------------------------------------------------------------------------
 * File:    autoreg.h
 * Date:    27/07/2011
 * Author:  Oscar Esteban oesteban@die.upm.es
 * Version: 0.1
 * License: BSD
 * --------------------------------------------------------------------------------------

 Copyright (c) 2011, Oscar Esteban - oesteban@die.upm.es
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


#ifndef AUTOREG_H_
#define AUTOREG_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>


#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientImageFilter.h>

#define USE_ADAPTATIVE_BINING
#include "MIRegistration.h"

namespace bpo = boost::program_options;

typedef float PixelRange;
const unsigned int Dimension = 3;
typedef itk::Image<PixelRange, Dimension> FixedImage;
typedef itk::Image<PixelRange, Dimension> MovingImage;

const double GAUSSIAN_FILTER_DEVIATION    =   4.0;
const unsigned int NUM_SAMPLES            = 50000;
const unsigned int MRES_LEVELS            =     2;
const unsigned int MI_BINS                =    70;


const unsigned int SPSA_MIN_ITERATIONS    =   100;
const unsigned int SPSA_NUM_PERTURBATIONS =     2;
const double       SPSA_INITIAL_STEP      =   1.0;

const double SPSA_A                       =  50.0; // 50
const double SPSA_C                       =   0.4;
const double SPSA_a                       = 1.0e4;
const unsigned int SPSA_MAX_ITERATIONS    =   800; // 500

const double GD_INITIAL_LEARNING_RATE     =  1000;
const unsigned int GD_MAX_ITERATIONS      =  1000;
const double GD_ALPHA                     =  0.95;
const double GD_a                         =  30.0;
const double GD_A                         =   100;


struct QuantileOption
{
	bool default_opt;
    boost::optional<float> lo_thres;
    boost::optional<float> hi_thres;
};

const QuantileOption FIXED_QUANTILES =  { true, 0.55, 0.999 }; // 0.55 - 0.75 // (0.65,0.97)
const QuantileOption MOVING_QUANTILES = { true, 0.60, 0.999 }; // (0.85, 0.999)

// Called by program_options to parse a set of Model arguments
void validate(boost::any& v, const std::vector<std::string>& values, QuantileOption*, int);

int main(int argc, char *argv[]);

#endif /* AUTOREG_H_ */
