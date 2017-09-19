/* sdsl - succinct data structures library
 *    Copyright (C) 2009 Simon Gog
 * 
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 * 
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see http://www.gnu.org/licenses/ .
 */
/*! \file rank_select_support_bp.hpp
 *    \brief TODO
 *    \author Tobias Heuer
 */
#ifndef INCLUDED_SDSL_RANK_SELECT_SUPPORT_BP
#define INCLUDED_SDSL_RANK_SELECT_SUPPORT_BP

#include <stack>
#include <limits>

#include "rank_support_v5.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{
    
    
    //! Class supports rank_1 and select_1 on BP sequences.
    /*!
     * \tparam t_sample_size If the depth of the BP sequence exceeds a certain threshold
     *              	     every t_sample_size-th position of a 1-bit is sampled.
     *  \tparam t_rank       Type of rank support used for the underlying BP sequence.
     *
     * \par Time complexity for select
     *        \f$ \Order{\log{\log{n}}} \f$
     * \par Space complexity:
     *        \f$ o(n) \f$. 
     * 
     */
    template<uint32_t t_sample_size = 1024>
    class rank_select_support_bp
    {
        
        const bit_vector*           m_v;
        rank_support_v5<>           m_bp_rank;
        int_vector<>                m_select_sample;
        bit_vector::value_type      m_max_excess;
        
        
        void copy(const rank_select_support_bp& rm) {
            m_v = rm.m_v;
            m_bp_rank = rm.m_bp_rank;
            m_bp_rank.set_vector(m_v);
            m_select_sample = rm.m_select_sample;
            m_max_excess = rm.m_max_excess;
        }
        
        
    private:
        
        //! Calculates the maximum depth of a leaf in the BP sequence.
        void calculate_maximum_excess_value() {
            m_max_excess = std::numeric_limits<bit_vector::value_type>::min();
            bit_vector::value_type cur_excess = 0;
            for(size_t i = 0; i < m_v->size(); ++i) {
                if((*m_v)[i]) cur_excess++;
                else cur_excess--;
                if(cur_excess > m_max_excess) m_max_excess = cur_excess;
            }
            if(m_v->size() <= 2) m_max_excess = 1;
        }        
        
        //! Samples the position of every t_sample_size-th 1 of the BP sequence. 
        void generate_select_sample() {
            size_t N = m_v->size()/2;
            m_select_sample = int_vector<>(N/t_sample_size+2,0);
            for(size_t i = 1; i < N; i += t_sample_size) {
                m_select_sample[i/t_sample_size] = select(i,0,false);
            }
            m_select_sample[N/t_sample_size+1] = select(N,0,false);
            util::bit_compress(m_select_sample);
        }
                                                                                 
                                                                                 
    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::value_type value_type;
        typedef rank_support_trait<1,1>  trait_type;
        
        //! Constructor
        /*
         *  \param v  Pointer to the bit vector of the BP sequence.
         */
        rank_select_support_bp(const bit_vector* v=nullptr) {
            if (v != nullptr) {
                set_vector(v);
                m_bp_rank = rank_support_v5<>(v);
                calculate_maximum_excess_value();
                //Construct sample datastructure, if depth of BP sequence is to big
                if(m_max_excess > t_sample_size) {
                    generate_select_sample();
                }
            }
        }
        
        //! Copy constructor
        rank_select_support_bp(const rank_select_support_bp& rm) {
            *this = rm;
        }
        
        //! Move constructor
        rank_select_support_bp(rank_select_support_bp&& rm) {
            *this = std::move(rm);
        }
        
        rank_select_support_bp& operator=(const rank_select_support_bp& rm) {
            if (this != &rm) {
                m_v = rm.m_v;
                m_bp_rank = rm.m_bp_rank;
                m_bp_rank.set_vector(m_v);
                m_select_sample = rm.m_select_sample;
                m_max_excess = rm.m_max_excess;
            }
            return *this;
        }
        
        rank_select_support_bp& operator=(rank_select_support_bp&& rm) {
            if (this != &rm) {
                m_v = std::move(rm.m_v);
                m_bp_rank = std::move(rm.m_bp_rank);
                m_bp_rank.set_vector(m_v);
                m_select_sample = std::move(rm.m_select_sample);
                m_max_excess = std::move(rm.m_max_excess);
            }
            return *this;
        }
        
        /*! Calculates the excess value at index i.
         * \param i The index of which the excess value should be calculated.
         */
        inline size_type excess(size_type idx) const {
            return (m_bp_rank(idx+1)<<1)-(idx+1);
        }
        
        /*! Returns the number of opening parentheses up to and including index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
         */
        inline size_type rank(size_type idx) const {
            return m_bp_rank.rank(idx);
        }
        
        /*! Returns the index of the i-th opening parenthesis.
         * \param i Number of the parenthesis to select.
         * \pre{ \f$1\leq i < rank(size())\f$ }
         * \post{ \f$ 0\leq select(i) < size() \f$ }
         */
        inline size_type select(size_type idx, size_type left = 0, bool use_select_samples=true) const {
            
            //Initial restriction interval [l,r] for select_1(idx)
            size_type l = std::max(left, 2*(idx-1)-std::min(static_cast<size_type>(m_max_excess),static_cast<size_type>(2*(idx-1))));
            size_type r = 2*(idx-1); 
            
            
            //If the Depth of the BP-Sequence is greater than t_sample_size, we use select use_select_samples
            //to further restrict the interval and speed up the running time
            if(use_select_samples && m_max_excess > t_sample_size) {
                size_t sample_idx = (idx-1)/t_sample_size;
                l = std::max(l,m_select_sample[sample_idx]);
                r = std::min(r,m_select_sample[sample_idx+1]);
            }
            
            size_type select_pos = 0, cur_rank = 0;
            
            //If the restriction interval of the i-th opening parenthesis is greater than 128,
            //we use the internal structure of rank_support_v5 to further restrict this interval
            //to size smaller than 64.
            if(r-l >= 128) {
                size_type bp_size = m_v->size();
                size_type N = m_bp_rank.m_basic_block.size();
                size_type start_pos = ((l>>10)&0xFFFFFFFFFFFFFFFEULL);
                select_pos = (start_pos>>1)<<11;
                
                //Prefix Sum of the 2048 bit superblock which contains select_1(idx)
                const uint64_t* p = m_bp_rank.m_basic_block.data() + start_pos;// (idx/2048)*2
                while(start_pos + 2 < N && *(p + 2) < idx) { p += 2; start_pos += 2; select_pos += 2048; }
                cur_rank = *p;
                
                //Calculate the prefix sum of the 348 bit basic block, which contains select_1(idx)
                //First  Block from bit 48 to bit 59
                //Second Block from bit 36 to bit 47
                //Third  Block from bit 24 to bit 35
                //Fourth Block from bit 12 to bit 23
                //Fifth  Block from bit  0 to bit 11
                const uint64_t basic_block_sum = *(p+1);
                size_type i = 1;
                for(; i < 6; ++i) {
                    size_type tmp_cur_rank = cur_rank + ((basic_block_sum>>(60-12*i))&0x7FFULL);
                    if(tmp_cur_rank >= idx || select_pos + ((i*6)<<6) >= bp_size) {
                        break;
                    }
                }
                select_pos += ((i-1)*6)<<6;
                cur_rank += (basic_block_sum>>(60-12*(i-1)))&0x7FFULL;
                
                //Calculate the prefix sum of the 64 bits block, which contains select_1(idx)
                const uint64_t* data = m_v->data()+(select_pos>>6);
                for(i = 0; i < 6; ++i) {
                    cur_rank += trait_type::full_word_rank(data,(i<<6)); 
                    if(cur_rank >= idx) {
                        break;
                    }
                }
                select_pos += i*64;
                cur_rank -=  trait_type::full_word_rank(data,(i<<6));
            } 
            //If the restriction interval of the i-th opening parenthesis is smaller than 128
            //and bigger than 64, we use simple binary search with the rank datastructure to 
            //further restrict the interval.
            else if(r-l >= 64) {
                while(r - l >= 64) {
                    size_type m = (l+r)/2;
                    if(m_bp_rank.rank(m) < idx) l = m+1;
                    else r = m;
                }
                select_pos = l;
                cur_rank = m_bp_rank(select_pos);
            }
            //If the restriction interval of the i-th opening parenthesis is smaller than 64,
            //we simple calculate the rank on position l.
            else {
                select_pos = l;
                cur_rank = m_bp_rank(select_pos);
            }
            
            //Calculate the exact position of select_1(idx) with binary search in a 64-bit word
            l = 0; r = 64;
            value_type data = m_v->get_int(select_pos,r-l);
            while(l < r) {
                size_type m = (l+r)/2;
                if(cur_rank + bits::cnt(data & bits::lo_set[m]) < idx) l = m+1;
                else r = m;
            }
            select_pos += (l - 1);
            
            return select_pos;
        }
        
        
        void swap(rank_select_support_bp& rm) {
            m_bp_rank.swap(rm.m_bp_rank); 
            m_select_sample.swap(rm.m_select_sample);
            m_max_excess = rm.m_max_excess;
        }
        
        
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_bp_rank.serialize(out, child, "bp_rank");
            written_bytes += m_select_sample.serialize(out, child, "select_sample");
            written_bytes += write_member(m_max_excess, out, child, "max_depth");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }
        
        void load(std::istream& in, const bit_vector* v=nullptr) {
            m_bp_rank.load(in,v);
            m_select_sample.load(in);
            read_member(m_max_excess,in);
            set_vector(v);
        }
        
        void set_vector(const bit_vector* v=nullptr) {
            m_v = v;
            m_bp_rank.set_vector(m_v);
        }
        
    };
    
} // end namespace sdsl
#endif
