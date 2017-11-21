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
/*! \file rmq_succinct_rec_new.hpp
 *    \brief rmq_succinct_rec_new.hpp contains the class rmq_succinct_rec_new which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
 *    \author Tobias Heuer
 */
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_REC_NEW
#define INCLUDED_SDSL_RMQ_SUCCINCT_REC_NEW

#include <stack>
#include <limits>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bits.hpp"
#include "bp_support_sada.hpp"
#include "bp_support_algorithm.hpp"
#include "rank_select_support_bp.hpp"
#include "suffix_tree_helper.hpp"
#include "util.hpp"


//! Namespace for the succinct data structure library.
namespace sdsl
{

template<bool t_min = true, uint32_t t_super_block_size=1024, uint32_t... t_block_sizes>
class rmq_succinct_rec_new;

template<bool t_min = false, uint32_t t_super_block_size=1024, uint32_t... t_block_sizes>
struct range_maximum_rec_new {
    typedef rmq_succinct_rec_new<t_min, t_super_block_size, t_block_sizes...> type;
};


//! A class to support range minimum or range maximum queries on a random access container.
/*!
 *  \tparam t_min        Specifies whether the data structure should answer range min/max queries (mimumum=true)
 *  \tparam t_bp_support Type of Support structure for the BPS-SCT.
 *
 * \par Time complexity
 *        \f$ \Order{1} \f$ for the range minimum/maximum queries if the balanced parentheses support structure supports constant time operations.
 * \par Space complexity:
 *        \f$ \Order{2n}+o(n) \f$ bits for the data structure ( \f$ n=size() \f$ ).
 *
 * \par Reference
 * H . Ferrada and G. Navarro.
 * Improved Range Minimum Queries.
 * In Proceedings of Data Compression Conference, DCC'16.
 *
 */
template<bool t_min, uint32_t t_super_block_size, uint32_t... t_block_sizes>
class rmq_succinct_rec_new
{
        using recursive_rmq = rmq_succinct_rec_new<t_min, t_block_sizes...>;
        using sparse_table = rmq_support_sparse_table<true,false>;

        bool                        m_use_sparse_rmq;       // Indicate, if sparse rmq derminates the recursion
        bit_vector                  m_gct_bp;               // BP-Sequence of the cartesian tree
        int_vector<>                m_min_excess_idx;       // Relative indicies of the minimum excess values inside a block
        int_vector<>                m_min_excess;           // Minimum excess values inside a block
        sparse_table                m_sparse_rmq;           // Sparse-RMQ (only build, if we derminate the recursion)
        recursive_rmq*              m_rmq_recursive;        // Recursive-RMQ
        rank_select_support_bp<>    m_rank_select;          // Rank and Select-Datastructure
        bit_vector::value_type      m_max_excess_v;         // Depth of the Cartesian Tree build over the original array
        bit_vector::value_type      m_max_excess_reverse_v; // Depth of the Cartesian Tree build over the reverse array

        void copy(const rmq_succinct_rec_new& rm) {
            m_use_sparse_rmq = rm.m_use_sparse_rmq;
            if (!m_use_sparse_rmq) {            
                m_gct_bp = rm.m_gct_bp;
                m_rank_select = rm.m_rank_select;
                m_rank_select.set_vector(&m_gct_bp);
                m_min_excess = rm.m_min_excess;
                m_min_excess_idx = rm.min_excess_idx;
                m_rmq_recursive = rm.m_rmq_recursive;
                m_max_excess_v = rm.m_max_excess_v;
                m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
            } else {
                m_min_excess = rm.m_min_excess;
                m_sparse_rmq = rm.m_sparse_rmq;
                m_sparse_rmq.set_vector(&m_min_excess);
            }
        }

    private:

        template<class t_rac>
        rmq_succinct_rec_new() : m_rmq_recursive(nullptr) { }

        /*!
         * We support several generalized cartesian tree layouts to minimize the depth of the tree.
         * All generalized cartesian tree, constructed with this method, use the leftmost path mapping. 
         * However, if we reverse the input sequence, we can build cartesian trees, which are isomorphic 
         * to cartesian trees with the rightmost path mapping. This has the advantage that we can dynamically
         * adapt the depth of cartesian tree during runtime to speed up depth dependent operations.         
         */
        /*!
         * \tparam t_strict  If true, the cartesian tree returns the leftmost minimum of sequence v
         *                   and the rightmost otherwise.
         * \tparam t_reverse If true, the cartesian tree is constructed over the reverse sequence of v.
         * \tparam t_rac     Random Access container
         *
         * \param v Sequence over which the cartesian tree should be constructed.
         * \param write_bp_sequence If true, the BP-Sequence of the cartesian tree is stored in
         *                          m_gct_bp. One can set this parameter to false, in order to
         *                          calculate only the depth of the cartesian tree.
         * \return The depth of the constructed cartesian tree.
         */
        template<bool t_strict, bool t_reverse, class t_rac>
        bit_vector::value_type construct_generalized_cartesian_tree(const t_rac* v, bool write_bp_sequence=true) {
            typedef min_max_trait<t_rac, true, t_strict> mm_trait;
            bit_vector::value_type max_excess = 0, cur_excess = 0;
            if (v->size() > 0) {
                long int cur_pos = t_reverse ? v->size()-1 : 0;
                size_t bp_cur_pos = 0;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min());
                if (write_bp_sequence) m_gct_bp[bp_cur_pos++] = 1;
                while ((!t_reverse && cur_pos <  ((long int) v->size())) || (t_reverse && cur_pos >= 0)) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos];
                    cur_pos += t_reverse ? -1 : 1;
                    while (mm_trait::compare(cur_elem, s.top()) && s.size() > 1) {
                        s.pop();
                        bp_cur_pos++; cur_excess--;
                    }
                    if (write_bp_sequence) m_gct_bp[bp_cur_pos++] = 1;
                    cur_excess++;
                    if (cur_excess > max_excess) max_excess = cur_excess;
                    s.push(cur_elem);
                }
                while (!s.empty()) {
                    s.pop();
                    bp_cur_pos++;
                }
            }
            return max_excess;
        }

        /*!
         * Constructs the recursive datastructure for our RMQ. The recursive RMQ is build over
         * the minimum excess values of the BP-Sequence of blocks of size t_super_block_size. 
         */
        void build_rmq_recursive() {
            size_type bp_size = m_gct_bp.size();
            size_type excess_block_size = bp_size/t_super_block_size + (bp_size % t_super_block_size != 0 ? 1 : 0); 
            m_min_excess = int_vector<>(excess_block_size,0);
            m_min_excess_idx = int_vector<>(excess_block_size,0);
            bit_vector::difference_type min_rel_ex = 0;
            for (size_t i = 0; i*t_super_block_size < bp_size; ++i) {
                uint64_t min_idx = near_rmq(m_gct_bp,i*t_super_block_size, std::min((i+1)*t_super_block_size - 1,bp_size-1),min_rel_ex);
                m_min_excess_idx[i] = min_idx-i*t_super_block_size;
                //Note: Minimum excess values are stored in reverse order, because our RMQ returns
                //      always the leftmost minimum, but in recursion we need the rightmost minimum.
                //      Reversing the sequence and reversing again the result insure that the rightmost
                //      minimum is returned.
                m_min_excess[excess_block_size - 1 - i] = m_rank_select.excess(min_idx);
            }
            m_rmq_recursive = new recursive_rmq(&m_min_excess);
            m_use_sparse_rmq = false;
            util::bit_compress(m_min_excess);
            util::bit_compress(m_min_excess_idx);
        }

        //! Determine rightmost minimum of three excess values
        inline std::pair<bit_vector::size_type, int_vector<>::value_type> 
               rightmost_minimum(const bit_vector::size_type    min_left_excess_idx, 
                                 const bit_vector::size_type    min_block_excess_idx, 
                                 const bit_vector::size_type    min_right_excess_idx,
                                 const int_vector<>::value_type min_left_excess, 
                                 const int_vector<>::value_type min_block_excess, 
                                 const int_vector<>::value_type min_right_excess) const {
            assert(min_left_excess_idx <= min_block_excess_idx); 
            assert(min_block_excess_idx <= min_right_excess_idx);
            if (min_right_excess <= min_left_excess && min_right_excess <= min_block_excess) {
                return std::make_pair(min_right_excess_idx,min_right_excess);
            } else if(min_block_excess <= min_left_excess) {
                return std::make_pair(min_block_excess_idx,min_block_excess);
            } else {
                return std::make_pair(min_left_excess_idx,min_left_excess);
            } 
        }

        /*!
         * The index of the minimum excess value of a block is stored relative 
         * to the block beginning. Therefore the real index for a block i is
         * calculated by adding i * t_super_block_size to m_min_excess_idx[i].
         */
        /*!
         * \param i Requested block for minimum excess index query.
         * \return The index of the minimum excess value in block i.
         */
        inline int_vector<>::value_type get_min_excess_idx(const size_t i) const {
            return m_min_excess_idx[i] + i * t_super_block_size;
        }

        //! Maps a block to its corresponding index in the minimum excess array (is stored in reverse order).
        /*!
         * \param i Block which should be mapped.
         * \return Index of i in the reversed minimum excess array.
         */
        inline int_vector<>::size_type map_to_min_excess(const size_t i) const {
            return m_min_excess.size()-1-i;
        }

        //! Returns the minimum excess value of a block i.
        /*!
         * \param i Requested block for minimum excess query. 
         * \return The minimum excess value of the requested block i.
         */
        inline int_vector<>::size_type get_min_excess(const size_t i) const {
            return m_min_excess[map_to_min_excess(i)];
        }

        /*!
         * The cartesian tree is either built on the original sequence or on the reverse.
         * Based on the type of the cartesian tree an index i has to be mapped back to the
         * original sequence (if reverse sequence is used).
         */
        /*!
         * \param i Index in the sequence over which the cartesian tree has been built.
         * \return Index in the original sequence.
         */
        inline int_vector<>::size_type map_index(const size_t i) {
            if(m_max_excess_v < m_max_excess_reverse_v) {
                return i;
            } else {
                return size() - 1 - i;
            }
        }

        /*!
         * Scan method to determine the rightmost minimum excess index
         * in the BP-Sequence. Is used to avoid second select operation, if
         * the corresponding opening parenthesis of an array index r is supposed
         * in 64-bit word, which start at position select_1_l.
         */
        /*!
         * \param l           Left index of the query.
         * \param r           Right index of the query.
         * \param select_1_l  Position of the opening parenthesis of l in the BP-Sequence.
         * \param min_idx     Stores the index of the minimum excess value in the BP-Sequence.
         * \return            True, if minimum excess value is found in a 64-bit word.
         */
        inline bool fast_rmq_scan(const bit_vector::size_type  l, 
                                  const bit_vector::size_type  r, 
                                  const bit_vector::size_type  select_1_l,
                                  int_vector<>::size_type& min_idx) const {
            value_type data = m_gct_bp.get_int(select_1_l+1);
            uint64_t one_cnt = bits::cnt(data)-1;
            //Case where the opening parenthesis of r is in data.
            //=> Simply 64-bit word scan determines the index with 
            //   minimum excess value in the BP-Sequence.
            if (l + one_cnt >= r) {
                size_t cur_select = l+1; min_idx = select_1_l;
                int cur_excess = ((l+1) << 1)-(select_1_l+1), min_excess = cur_excess;
                for (size_t k = select_1_l+1; k <= select_1_l+65; ++k) {
                    if (m_gct_bp[k]) {
                        cur_excess++; cur_select++;
                        if (cur_select == r + 2) break;
                    } else {
                        cur_excess--;
                        if (cur_excess <= min_excess) {
                            min_excess = cur_excess;
                            min_idx = k;
                        }
                    }
                }
                min_idx = (min_excess+min_idx)>>1;
                return true;
            }
            else {
                return false;
            }
        }

    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::value_type value_type;
        typedef typename int_vector<>::value_type i_value_type;

        //! Default constructor
        rmq_succinct_rec_new() {}

        //! Construct the recursive RMQ datastructure.
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_rec_new(const t_rac* v=nullptr) : m_rmq_recursive(nullptr) {
            if (v != nullptr) {
                size_type bp_size = 2*v->size()+2;
                m_gct_bp = bit_vector(bp_size,0);
                if(t_super_block_size > 0 && t_super_block_size <= bp_size) {
                    m_max_excess_v = construct_generalized_cartesian_tree<true,false>(v,false);
                    m_max_excess_reverse_v = construct_generalized_cartesian_tree<false,true>(v,false);
                    //Choose Cartesian Tree with minimal depth.
                    if (m_max_excess_v < m_max_excess_reverse_v) {
                        construct_generalized_cartesian_tree<true,false>(v);
                    } else {
                        construct_generalized_cartesian_tree<false,true>(v);
                    }
                    m_rank_select = rank_select_support_bp<>(&m_gct_bp);
                    build_rmq_recursive();
                }
                else {
                    //In case of building the Sparse-RMQ, we need to store the input array
                    //to serialize the datastructure.
                    m_min_excess = int_vector<>(v->size(),0);
                    for (size_t i = 0; i < v->size(); ++i) 
                        m_min_excess[i] = (*v)[i];
                    m_sparse_rmq = sparse_table(&m_min_excess);
                    m_use_sparse_rmq = true;
                }
            }
        }

        //! Destructor
        ~rmq_succinct_rec_new() {
            if (t_super_block_size > 0 && m_gct_bp.size() > 2) {
                delete m_rmq_recursive;
            }
        }

        //! Copy constructor
        rmq_succinct_rec_new(const rmq_succinct_rec_new& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_rec_new(rmq_succinct_rec_new&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_rec_new& operator=(const rmq_succinct_rec_new& rm) {
            if (this != &rm) {
                m_use_sparse_rmq = rm.m_use_sparse_rmq;
                if (!m_use_sparse_rmq) {
                    m_gct_bp = rm.m_gct_bp; 
                    m_rank_select = rm.m_rank_select;
                    m_rank_select.set_vector(&m_gct_bp);
                    m_min_excess = rm.m_min_excess;
                    m_min_excess_idx = rm.min_excess_idx;
                    m_rmq_recursive = rm.m_rmq_recursive;
                    m_max_excess_v = rm.m_max_excess_v;
                    m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
                } else {
                    m_min_excess = rm.m_min_excess;
                    m_sparse_rmq = rm.m_sparse_rmq;
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        rmq_succinct_rec_new& operator=(rmq_succinct_rec_new&& rm) {
            if (this != &rm) {
                m_use_sparse_rmq = rm.m_use_sparse_rmq;
                if (!m_use_sparse_rmq) {
                    m_gct_bp = std::move(rm.m_gct_bp);
                    m_rank_select = std::move(rm.m_rank_select);
                    m_rank_select.set_vector(&m_gct_bp);
                    m_min_excess = std::move(rm.m_min_excess);
                    m_min_excess_idx = std::move(rm.m_min_excess_idx);
                    m_rmq_recursive = std::move(rm.m_rmq_recursive);
                    m_max_excess_v = rm.m_max_excess_v;
                    m_max_excess_reverse_v = rm.m_max_excess_reverse_v;
                } else {
                    m_min_excess = std::move(rm.m_min_excess);
                    m_sparse_rmq = std::move(rm.m_sparse_rmq);
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        void swap(rmq_succinct_rec_new& rm) {
            std::swap(m_use_sparse_rmq, rm.m_use_sparse_rmq);
            if (!m_use_sparse_rmq) {
                m_gct_bp.swap(rm.m_gct_bp);
                util::swap_support(m_rank_select, rm.m_rank_select,
                                &m_gct_bp, &(rm.m_gct_bp));
                m_min_excess.swap(rm.m_min_excess);
                m_min_excess_idx.swap(rm.m_min_excess_idx);
                *m_rmq_recursive.swap(*rm.m_rmq_recursive);
                std::swap(m_max_excess_v, rm.m_max_excess_v);
                std::swap(m_max_excess_reverse_v, rm.m_max_excess_reverse_v);
            } else {
                m_min_excess.swap(rm.m_min_excess);
                util::swap_support(m_sparse_rmq, rm.m_sparse_rmq,
                                   &m_min_excess, &(rm.min_excess));
            }
        }


        //! Range minimum/maximum query for the supported random access container v.
        /*!
         * \param l Leftmost position of the interval \f$[\ell..r]\f$.
         * \param r Rightmost position of the interval \f$[\ell..r]\f$.
         * \return The minimal index i with \f$\ell \leq i \leq r\f$ for which \f$ v[i] \f$ is minimal/maximal.
         * \pre
         *   - r < size()
         *   - \f$ \ell \leq r \f$
         * \par Time complexity
         *      \f$ \Order{1} \f$
         */
        size_type operator()(const size_type l, const size_type r) {
            assert(l <= r); assert(r < size());
            if (l == r) return l;
            if(m_use_sparse_rmq) {
                 return m_sparse_rmq(l,r);
            }

            size_type tmp_l = map_index(l), tmp_r = map_index(r);
            if (tmp_l > tmp_r) std::swap(tmp_l,tmp_r);

            size_type i = m_rank_select.select(tmp_l+2)-1;

            //If the interval is small, we can try to scan the interval to avoid second select.
            if (tmp_r - tmp_l < 64) {
                size_type min_idx = 0;
                if(fast_rmq_scan(tmp_l,tmp_r,i,min_idx)) {
                    return map_index(min_idx);
                }
            }
            size_type j = m_rank_select.select(tmp_r+2,i);

            size_type block_i = (i+t_super_block_size-1)/t_super_block_size;
            size_type block_j = j/t_super_block_size;
            
            bit_vector::difference_type min_rel_ex = 0;
            //Case: Query spans less than three blocks
            if (block_i >= block_j) {
                block_i = block_i-(i != 0);
                size_type min_excess_idx = i;
                if (block_i == block_j) {
                    min_excess_idx = get_min_excess_idx(block_i);
                } else {
                    size_type min_left_excess_idx = get_min_excess_idx(block_i);
                    size_type min_right_excess_idx = get_min_excess_idx(block_j);
                    min_excess_idx = (get_min_excess(block_i) < get_min_excess(block_j) ? min_left_excess_idx : min_right_excess_idx);
                }
                if (min_excess_idx < i || min_excess_idx > j) {
                    min_excess_idx = near_rmq(m_gct_bp,i,j-1,min_rel_ex); 
                }
                return map_index(m_rank_select.rank(min_excess_idx+1)-1);
            } 
            //Case: Query can be divided into three queries [l,l'), [l',r'] and (r',r], where [l',r'] is
            //      aligned with the blocks. 
            else {
                //Block aligned query [l',r']
                size_type min_block_idx = (*m_rmq_recursive)(map_to_min_excess(block_j-1), map_to_min_excess(block_i));
                //Note: The recursive RMQ is built over the reverse minimum excess array and returns the 
                //      leftmost minimum (which corresponds to rightmost into the forward array). Accessing 
                //      the minimum excess indicies requires to map min_block_idx to the forward minimum excess
                //      array with map_to_min_excess. 
                size_type min_block_excess_idx = get_min_excess_idx(map_to_min_excess(min_block_idx));
                i_value_type min_block_excess = m_min_excess[min_block_idx];
                
                //Query [l,l') for the left block 
                size_type min_left_excess_idx = get_min_excess_idx(block_i-(i != 0)); 
                i_value_type min_left_excess = get_min_excess(block_i-(i != 0));
                if (min_left_excess < min_block_excess && min_left_excess_idx < i) {
                    min_left_excess_idx = near_rmq(m_gct_bp,i,t_super_block_size*block_i,min_rel_ex); 
                    min_left_excess = m_rank_select.excess(min_left_excess_idx); 
                }

                //Query (r',r] for the right block
                size_type min_right_excess_idx = get_min_excess_idx(block_j); 
                i_value_type min_right_excess = get_min_excess(block_j);
                if (min_right_excess <= min_block_excess && min_right_excess_idx > j) {
                    min_right_excess_idx = near_rmq(m_gct_bp,t_super_block_size*block_j,j,min_rel_ex); 
                    min_right_excess = m_rank_select.excess(min_right_excess_idx);
                }

                auto rmq_min = rightmost_minimum(min_left_excess_idx,min_block_excess_idx,min_right_excess_idx,
                                                 min_left_excess,min_block_excess,min_right_excess);
                size_type min_idx = rmq_min.first;
                int_vector<>::value_type min_ex = rmq_min.second;
                return map_index((min_ex+min_idx)>>1);
            }
        }


        size_type size()const {
            if(m_use_sparse_rmq) {
                return m_sparse_rmq.size();
            }
            else {
              return (m_gct_bp.size()-2)/2;
            }
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_use_sparse_rmq, out, child, "m_use_sparse_rmq");
            if(m_rmq_recursive != nullptr) {
                written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
                written_bytes += m_rank_select.serialize(out, child, "rank_select_bp");
                written_bytes += m_min_excess.serialize(out, child, "min_excess");
                written_bytes += m_min_excess_idx.serialize(out, child, "min_excess_idx");
                written_bytes += m_rmq_recursive->serialize(out, child, "rmq_recursive");
                written_bytes += write_member(m_max_excess_v, out, child, "m_max_excess_v");
                written_bytes += write_member(m_max_excess_reverse_v, out, child, "m_max_excess_v_reverse");
            } else {
                written_bytes += m_min_excess.serialize(out, child, "min_excess");
                written_bytes += m_sparse_rmq.serialize(out, child, "sparse_rmq");
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            read_member(m_use_sparse_rmq,in);
            if (!m_use_sparse_rmq) {
                m_gct_bp.load(in);
                m_rank_select.load(in, &m_gct_bp);
                m_min_excess.load(in);
                m_min_excess_idx.load(in);
                m_rmq_recursive = new recursive_rmq();
                m_rmq_recursive->load(in);
                read_member(m_max_excess_v,in);
                read_member(m_max_excess_reverse_v,in);
            } else {
                m_min_excess.load(in);
                m_sparse_rmq.load(in,&m_min_excess);
            }
        }
};

} // end namespace sdsl
#endif

