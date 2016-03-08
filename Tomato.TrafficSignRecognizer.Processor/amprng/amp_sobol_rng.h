/*----------------------------------------------------------------------------
 * Copyright (c) Microsoft Corp. 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not 
 * use this file except in compliance with the License.  You may obtain a copy 
 * of the License at http://www.apache.org/licenses/LICENSE-2.0  
 *
 * THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY 
 * KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED 
 * WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE, 
 * MERCHANTABLITY OR NON-INFRINGEMENT. 
 *
 * See the Apache Version 2.0 License for specific language governing 
 * permissions and limitations under the License.
 *  
 * The data set of direction numbers for Sobol sequences is based the project 
 * on or incorporating material from the following project(s):
 * Sobol sequence generator, available at http://web.maths.unsw.edu.au/~fkuo/sobol/index.html
 * 
 * Specifically it is based on the file "new-joe-kuo-6.21201" containing the direction 
 * numbers obtained using the search criterion D(6) up to dimension 21201.
 *
 * Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 *  * Neither the names of the copyright holders nor the names of the
 *    University of New South Wales and the University of Waikato
 *    and its contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * -----------------------------------------------------------------------------------
 * 
 * File: amp_sobol_rng.h
 * 
 * Implements 32b Sobol quasi random number generator using C++ AMP
 *------------------------------------------------------------------------------------ */

#pragma once
#ifndef _AMP_SOBOL_RNG_H
#define _AMP_SOBOL_RNG_H

#include "amp_rand_collection.h"
#include "xxamp_sobol_rng.h"

/// This is the class implementing quasi_sobol sequence engine
template<unsigned sobol_dimension>
class sobol_rng
{
    static_assert(sobol_dimension >= 1 && 
                  sobol_dimension <= sobol_rng_lib::dimension_limit, 
                  "the specified sobol RNG dimension is not supported");
    
    typedef sobol_rng_lib::direction_num_view direction_num_view;

public:

    // Multi-dimension sobol number
    template<typename value_type>
    class sobol_number
    {
    private:
        value_type values[sobol_dimension];

    public:
        const value_type& operator[](int dim) const restrict(cpu, amp)
        {
            return values[dim];
        }

        value_type& operator[](int dim) restrict(cpu, amp)
        {
            return values[dim];
        }
    };

    /// Setup state
    void initialize(const direction_num_view &direction_nums, unsigned skipahead = 0) restrict(cpu, amp)
    {
        state.position = 0;

        for (int i=0; i<sobol_dimension; i++)
        {
            state.sobol_num[i] = 0;
        }

        skip(direction_nums, skipahead);
    }

    // Skip ahead from current position
    void skip(const direction_num_view &direction_nums, unsigned n) restrict(cpu, amp)
    {
        if (n == 0) return;

        // Update position. Overflow is handled as position wrapped to zero. 
        state.position += n;

        // Compute gray code from position
        unsigned gray_code = state.position ^ (state.position>>1);
        
        // Compute new state 
        for (int j=0; j<sobol_rng_lib::rng_bits; ++j)
        {
            if (gray_code & (1<<j))
            {
                for (int i=0; i<sobol_dimension; ++i)
                {
                    state.sobol_num[i] ^= direction_nums(i,j);
                }
            }
        }
    }

    /// Get to the next state
    void next(const direction_num_view &direction_nums) restrict(cpu, amp)
    {
        // Update position. Overflow is handled as position wrapped to zero. 
        state.position++;

        // Find the rightmost zero digit
        unsigned j = 0;
        unsigned p = state.position;
        while (p & 1) 
        {
            p >>= 1;
            j++;
        }

        // Compute new state 
        for (int i=0; i<sobol_dimension; ++i)
        {
            state.sobol_num[i] ^= direction_nums(i,j);
        }
    }

    /// Get the unsigned sobol number at specified dimension
    unsigned get_uint(unsigned dim) restrict(cpu, amp)
    {
        return state.sobol_num[dim-1];
    }

    /// Get the floating point sobol number between [0, 1) at specified dimension
    float get_single(unsigned dim) restrict(cpu, amp)
    {
        // Trick to avoid overflow of (1<<32)
        const float scale_factor = 0.5f/(unsigned)(1<<(sobol_rng_lib::rng_bits-1)); 

        return static_cast<float>(get_uint(dim))* scale_factor;
    }

    /// Get the floating point random number between [1, 2) at specified dimension
    float get_single12(unsigned dim) restrict(cpu, amp)
    {
        return get_single(dim) + 1.0f;
    }

private:
    struct sobol_state
    {
        unsigned               position;
        sobol_number<unsigned> sobol_num;
    };

    sobol_state state;
};

/// This class initialize all the Sobol sequence RNG engines
template<typename rng_type, int rank>
class  sobol_rng_collection : public amp_rand_collection<rng_type, rank>
{
private:
    // Direction number array_view
    concurrency::array_view<unsigned, 2> direction_num_av;

    sobol_rng_collection()
    {
    }

public:
    sobol_rng_collection(const concurrency::extent<rank> rand_extent, unsigned skipahead=0)
        : _base(rand_extent), 
          direction_num_av(concurrency::extent<2>(sobol_rng_lib::dimension_limit, sobol_rng_lib::rng_bits), 
                           sobol_rng_lib::direction_nums)
    {
        concurrency::array_view<rng_type, rank> rng_av(m_rng_av);
        concurrency::array_view<unsigned, 2> direction_nums(direction_num_av);

        parallel_for_each(rand_extent, [=] (concurrency::index<rank> idx) restrict(amp)
        {
            rng_av[idx].initialize(direction_nums, skipahead);
        });
    }

	sobol_rng_collection(const concurrency::accelerator_view& acc_view, const concurrency::extent<rank> rand_extent, unsigned skipahead = 0)
		:_base(acc_view, rand_extent),
		direction_num_av(concurrency::extent<2>(sobol_rng_lib::dimension_limit, sobol_rng_lib::rng_bits),
			sobol_rng_lib::direction_nums)
	{
		concurrency::array_view<rng_type, rank> rng_av(m_rng_av);
		concurrency::array_view<unsigned, 2> direction_nums(direction_num_av);

		parallel_for_each(acc_view, rand_extent, [=](concurrency::index<rank> idx) restrict(amp)
		{
			rng_av[idx].initialize(direction_nums, skipahead);
		});
	}

    // Get direction numbers
    const sobol_rng_lib::direction_num_view& direction_numbers() const restrict(cpu, amp)
    {
        return direction_num_av;
    }
};

#endif // _AMP_SOBOL_RNG_H