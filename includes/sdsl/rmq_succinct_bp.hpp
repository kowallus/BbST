/* sdsl - succinct data structures library
    Copyright (C) 2009 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file rmq_succinct_bp.hpp
    \brief rmq_succinct_bp.hpp contains the class rmq_succinct_bp which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
    \author Tobias Heuer
*/
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_BP
#define INCLUDED_SDSL_RMQ_SUCCINCT_BP

#include <stack>
#include <limits>

#include "rmq_support.hpp"
#include "int_vector.hpp"
#include "bp_support_sada.hpp"
#include "suffix_tree_helper.hpp"
#include "util.hpp"

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<bool t_min = true,
         class t_bp_support = bp_support_sada<256,32,rank_support_v5<> > >
class rmq_succinct_bp;

template<class t_bp_support = bp_support_sada<256,32,rank_support_v5<> > >
struct range_maximum_bp {
    typedef rmq_succinct_bp<false, t_bp_support> type;
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
template<bool t_min, class t_bp_support>
class rmq_succinct_bp
{
        bit_vector            m_gct_bp;         //!< A bit vector which contains the BP-GT of the input container.
        t_bp_support          m_gct_bp_support; //!< Support structure for the BPS-GT

        void copy(const rmq_succinct_bp& rm) {
            m_gct_bp = rm.m_gct_bp;
            m_gct_bp_support = rm.m_gct_bp_support;
            m_gct_bp_support.set_vector(&m_gct_bp);
        }

        
    private:
        
        template<class t_rac>
        void construct_generalized_cartesian_tree_rightmost(const t_rac* v) {
            if(v->size() > 0) {
                int64_t cur_pos = v->size()-1;
                int64_t bp_cur_pos = m_gct_bp.size()-1;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min()); 
                m_gct_bp[bp_cur_pos--] = 0;
                while(cur_pos >= 0) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos--];
                    while(s.top() >= cur_elem) {
                        s.pop();
                        m_gct_bp[bp_cur_pos--] = 1;
                    }
                    bp_cur_pos--;
                    s.push(cur_elem);
                }
                while(!s.empty()) {
                    s.pop();
                    m_gct_bp[bp_cur_pos--] = 1;
                }
            }
        }
        
        template<class t_rac>
        void construct_generalized_cartesian_tree_leftmost(const t_rac* v) {
            if(v->size() > 0) {
                size_t cur_pos = 0;
                size_t bp_cur_pos = 0;
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min()); 
                m_gct_bp[bp_cur_pos++] = 1;
                while(cur_pos < v->size()) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos++];
                    while(s.top() > cur_elem) {
                        s.pop();
                        bp_cur_pos++;
                    }
                    m_gct_bp[bp_cur_pos++] = 1;
                    s.push(cur_elem);
                }
                while(!s.empty()) {
                    s.pop();
                    bp_cur_pos++;
                }
            }
        }
        
        
    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::size_type value_type;
        typedef t_bp_support                   bp_support_type;

        const bit_vector&      gct_bp         = m_gct_bp;
        const bp_support_type& gct_bp_support = m_gct_bp_support;

        //! Default constructor
        rmq_succinct_bp() {}

        //! Constructor
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_bp(const t_rac* v=nullptr) {
            if (v != nullptr) {
                m_gct_bp = bit_vector(2*v->size()+2,0);
                construct_generalized_cartesian_tree_leftmost(v);
                m_gct_bp_support = bp_support_type(&m_gct_bp); 
            }
        }

        //! Copy constructor
        rmq_succinct_bp(const rmq_succinct_bp& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_bp(rmq_succinct_bp&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_bp& operator=(const rmq_succinct_bp& rm) {
            if (this != &rm) {
                m_gct_bp = rm.m_gct_bp;
                m_gct_bp_support = rm.m_gct_bp_support;
                m_gct_bp_support.set_vector(&m_gct_bp);
            }
            return *this;
        }

        rmq_succinct_bp& operator=(rmq_succinct_bp&& rm) {
            if (this != &rm) {
                m_gct_bp = std::move(rm.m_gct_bp);
                m_gct_bp_support = std::move(rm.m_gct_bp_support);
                m_gct_bp_support.set_vector(&m_gct_bp);              
            }
            return *this;
        }

        void swap(rmq_succinct_bp& rm) {
            m_gct_bp.swap(rm.m_gct_bp);
            util::swap_support(m_gct_bp_support, rm.m_gct_bp_support,
                               &m_gct_bp, &(rm.m_gct_bp));       
        }
        
        void print_bp() {
            size_t N = m_gct_bp.size();
            for(size_t i = 0; i < N; ++i) {
                std::cout << (m_gct_bp[i] ? '(' : ')');
            }
            std::cout << std::endl;
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
        size_type operator()(const size_type l, const size_type r)const {
            assert(l <= r); assert(r < size());
            size_type i     = m_gct_bp_support.select(l+2)-1;
            size_type j     = m_gct_bp_support.select(r+2);
            size_type rmq_e = m_gct_bp_support.rmq(i,j);
            return m_gct_bp_support.rank(rmq_e)-1;
        }

        size_type size()const {
            return (m_gct_bp.size()-2)/2;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
            written_bytes += m_gct_bp_support.serialize(out, child, "gct_bp_support");          
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_gct_bp.load(in);
            m_gct_bp_support.load(in, &m_gct_bp);
        }
};

} // end namespace sdsl
#endif
