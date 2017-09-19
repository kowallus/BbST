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
/*! \file rmq_succinct_rec_old.hpp
    \brief rmq_succinct_rec_old.hpp contains the class rmq_succinct_bp_fast which supports range minimum or range maximum queries on a random access container in constant time and \f$2 n+o(n) bits\f$ space.
    \author Tobias Heuer
*/
#ifndef INCLUDED_SDSL_RMQ_SUCCINCT_REC_OLD
#define INCLUDED_SDSL_RMQ_SUCCINCT_REC_OLD

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

// #define MEASURE_TIMINGS

#ifdef MEASURE_TIMINGS
using MeasurePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
MeasurePoint start, end;

#define NANO 1000000000
inline MeasurePoint measure_point()
{
    return std::chrono::high_resolution_clock::now();
}

double ns()
{
    std::chrono::duration<double> elapsed_seconds = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    return elapsed_seconds.count()*NANO;
}


class Stats
{
    private:
        using StatsMap = std::map<std::string, double>;

    public:
        Stats(const Stats&) = delete;
        Stats(Stats&&) = delete;
        Stats& operator= (const Stats&) = delete;
        Stats& operator= (Stats&&) = delete;

        void add(const std::string& key, double value) {
            _stats[key] += std::max((value-_avg_timer_overhead),0.0);
        }

        void addToTotal(const std::string& key, double value) {
            _stats[key] += std::max((value-(key != "Range" ? _avg_timer_overhead : 0.0)),0.0);
        }


        std::string toStringAndClear() {
            std::ostringstream s;
            for (auto& stat : _stats) {
                s << " " << stat.first << "=" << stat.second;
            }
            reset();
            return s.str();
        }

        void reset() {
            _stats.clear();
            for (std::string& key : keys) {
                addToTotal(key,0.0);
            }
        }

        void calculate_average_timer_overhead(size_t N) {
            double tmp_timer_overhead = 0.0;
            for (size_t i = 0; i < N; ++i) {
                start = measure_point();
                end = measure_point();
                tmp_timer_overhead += ns();
            }
            _avg_timer_overhead = tmp_timer_overhead/N;
        }

        static Stats& instance() {
            static Stats instance;
            return instance;
        }


    private:
        Stats() :
            _stats(), _avg_timer_overhead(0.0) { reset(); calculate_average_timer_overhead(100000); }
        StatsMap _stats;
        double _avg_timer_overhead;
        std::vector<std::string> keys = {"Rank","Sparse_RMQ","Scan","Select","min_excess","min_excess_idx","Other"};
};


#define TIME_MEASURE(X,key)                                \
    start = measure_point();                               \
    X                                                      \
    end = measure_point();                                 \
    Stats::instance().addToTotal(key,ns());                \
 
#define RETURN_TIME_MEASURE(type,X,key,t_block_size)    \
    start = measure_point();                                \
    type val = X                                            \
    end = measure_point();                                  \
    Stats::instance().addToTotal(key,ns());                 \
    if(t_block_size > 0)                                    \
        std::cout << "TIMING_RESULT" << Stats::instance().toStringAndClear() << std::endl; \
    return val;                                             \
 
#endif

#ifndef MEASURE_TIMINGS
#define TIME_MEASURE(X,key) X
#define RETURN_TIME_MEASURE(type,X,key,t_block_size) return X
#endif


//! Namespace for the succinct data structure library.
namespace sdsl
{

template<uint32_t t_super_block_size=512,uint32_t t_block_size=128, bool t_min = true>
class rmq_succinct_rec_old;

template<uint32_t t_super_block_size=512, uint32_t t_block_size=128>
struct range_maximum_bp_fast_rec_old {
    typedef rmq_succinct_rec_old<t_super_block_size, t_block_size, false> type;
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
template<uint32_t t_super_block_size, uint32_t t_block_size, bool t_min>
class rmq_succinct_rec_old
{
        bit_vector                  m_gct_bp;         //!< A bit vector which contains the BP-GT of the input container.
        int_vector<>                m_min_excess_idx;
        int_vector<>                m_min_excess;
        rmq_support_sparse_table<>  m_sparse_rmq;
        rmq_succinct_rec_old<t_block_size,0>*  m_rmq_recursive;
        rank_select_support_bp<>      m_rank_select;

        void copy(const rmq_succinct_rec_old& rm) {
            m_gct_bp = rm.m_gct_bp;
            m_rank_select = rm.m_rank_select;
            m_rank_select.set_vector(&m_gct_bp);
            m_min_excess = rm.m_min_excess;
            m_min_excess_idx = rm.min_excess_idx;
            if (t_block_size > 0) {
                m_rmq_recursive = rm.m_rmq_recursive;
            } else {
                m_sparse_rmq = rm.m_sparse_rmq;
                m_sparse_rmq.set_vector(&m_min_excess);
            }
        }


    private:
        template<class t_rac>
        rmq_succinct_rec_old() { }

        template<class t_rac>
        void construct_generalized_cartesian_tree_leftmost(const t_rac* v) {
            if (v->size() > 0) {
                size_t cur_pos = 0;
                size_t bp_cur_pos = 0;
                //TODO: Use sorted_multi_stack_support
                // Use positions instead of values in succinct_stack
                std::stack<typename t_rac::value_type> s;
                s.push(std::numeric_limits<typename t_rac::value_type>::min());
                m_gct_bp[bp_cur_pos++] = 1;
                while (cur_pos < v->size()) {
                    typename t_rac::value_type cur_elem = (*v)[cur_pos++];
                    if (t_block_size == 0) {
                        while (s.top() >= cur_elem && s.size() > 1) {
                            s.pop();
                            bp_cur_pos++;
                        }
                    } else {
                        while (s.top() > cur_elem && s.size() > 1) {
                            s.pop();
                            bp_cur_pos++;
                        }
                    }
                    m_gct_bp[bp_cur_pos++] = 1;
                    s.push(cur_elem);
                }
                while (!s.empty()) {
                    s.pop();
                    bp_cur_pos++;
                }
            }
        }

        void construct_minimum_excess_rmq() {
            size_type bp_size = m_gct_bp.size();
            m_min_excess = int_vector<>(bp_size/t_super_block_size+1,0);
            m_min_excess_idx = int_vector<>(bp_size/t_super_block_size+1,0);
            bit_vector::difference_type min_rel_ex = 0;
            for (size_t i = 0; i*t_super_block_size < bp_size; ++i) {
                uint64_t min_idx = near_rmq(m_gct_bp,i*t_super_block_size, std::min((i+1)*t_super_block_size - 1,bp_size-1),min_rel_ex);
                m_min_excess_idx[i] = min_idx-i*t_super_block_size;
                m_min_excess[i] = m_rank_select.excess(min_idx);
            }
            if (t_block_size > 0) {
                m_rmq_recursive = new rmq_succinct_rec_old<t_block_size,0>(&m_min_excess);
            } else {
                m_sparse_rmq = rmq_support_sparse_table<>(&m_min_excess);
            }
            util::bit_compress(m_min_excess);
            util::bit_compress(m_min_excess_idx);
        }

        inline std::pair<bit_vector::size_type, int_vector<>::value_type> min_ex(const bit_vector::size_type i1, const bit_vector::size_type i2, const bit_vector::size_type i3,
                const int_vector<>::value_type i1_ex, const int_vector<>::value_type i2_ex, const int_vector<>::value_type i3_ex) const {
            assert(i1 <= i2); assert(i2 <= i3);
            //std::cout << i1_ex << " " << i2_ex << " " << i3_ex << std::endl;
            if (i1_ex > i2_ex) {
                if (i2_ex < i3_ex) return std::make_pair(i2,i2_ex);
                else return  std::make_pair(i3,i3_ex);;
            } else if (i1_ex < i2_ex) {
                if (i1_ex < i3_ex) return  std::make_pair(i1,i1_ex);
                else return std::make_pair(i3,i3_ex);
            } else {
                if (i2_ex < i3_ex) return  std::make_pair(i2,i2_ex);
                else return  std::make_pair(i3,i3_ex);;
            }
        }

        inline int_vector<>::value_type get_min_excess_idx(size_t i) const {
            return m_min_excess_idx[i] + i * t_super_block_size;
        }





    public:
        typedef typename bit_vector::size_type size_type;
        typedef typename bit_vector::value_type value_type;
        typedef typename int_vector<>::value_type i_value_type;

        const bit_vector&                  gct_bp         = m_gct_bp;
        rmq_succinct_rec_old<t_block_size,0>* rmq_recursive = m_rmq_recursive;
        const rmq_support_sparse_table<>&  sparse_rmq     = m_sparse_rmq;
        const int_vector<>&                min_excess     = m_min_excess;
        const int_vector<>&                min_excess_idx = m_min_excess_idx;
        const rank_select_support_bp<>       rank_select    = m_rank_select;
        //const int_vector<>&                min_excess_idx64 = m_min_excess_idx64;


        size_t num_avoided_selects;
        size_t num_queries;


        //! Default constructor
        rmq_succinct_rec_old() {}

        //! Constructor
        /*! \tparam t_rac A random access container.
         *  \param  v     Pointer to container object.
         */
        template<class t_rac>
        rmq_succinct_rec_old(const t_rac* v=nullptr) {
            num_queries = 0;
            num_avoided_selects = 0;
            if (v != nullptr) {
                //TODO: Use construct_supercartesian_tree_bp_succinct in suffix_tree_helper.hpp
                m_gct_bp = bit_vector(2*v->size()+2,0);
                construct_generalized_cartesian_tree_leftmost(v);
                m_rank_select = rank_select_support_bp<>(&m_gct_bp);
                construct_minimum_excess_rmq();
            }
        }

        ~rmq_succinct_rec_old() {
            if (t_block_size > 0 && m_gct_bp.size() > 2) {
                delete m_rmq_recursive;
            }
        }

        //! Copy constructor
        rmq_succinct_rec_old(const rmq_succinct_rec_old& rm) {
            *this = rm;
        }

        //! Move constructor
        rmq_succinct_rec_old(rmq_succinct_rec_old&& rm) {
            *this = std::move(rm);
        }

        rmq_succinct_rec_old& operator=(const rmq_succinct_rec_old& rm) {
            if (this != &rm) {
                m_gct_bp = rm.m_gct_bp;
                m_rank_select = rm.m_rank_select;
                m_rank_select.set_vector(&m_gct_bp);
                m_min_excess = rm.m_min_excess;
                m_min_excess_idx = rm.min_excess_idx;
                if (t_block_size > 0) {
                    m_rmq_recursive = rm.m_rmq_recursive;
                } else {
                    m_sparse_rmq = rm.m_sparse_rmq;
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        rmq_succinct_rec_old& operator=(rmq_succinct_rec_old&& rm) {
            if (this != &rm) {
                m_gct_bp = std::move(rm.m_gct_bp);
                m_rank_select = std::move(rm.m_rank_select);
                m_rank_select.set_vector(&m_gct_bp);
                m_min_excess = std::move(rm.m_min_excess);
                m_min_excess_idx = std::move(rm.m_min_excess_idx);
                if (t_block_size > 0) {
                    m_rmq_recursive = std::move(rm.m_rmq_recursive);
                } else {
                    m_sparse_rmq = std::move(rm.m_sparse_rmq);
                    m_sparse_rmq.set_vector(&m_min_excess);
                }
            }
            return *this;
        }

        void swap(rmq_succinct_rec_old& rm) {
            m_gct_bp.swap(rm.m_gct_bp);
            util::swap_support(m_rank_select, rm.m_rank_select,
                               &m_gct_bp, &(rm.m_gct_bp));
            m_min_excess.swap(rm.m_min_excess);
            m_min_excess_idx.swap(rm.m_min_excess_idx);
            if (t_block_size > 0) {
                *m_rmq_recursive.swap(*rm.m_rmq_recursive);
            } else {
                m_sparse_rmq = std::move(rm.m_sparse_rmq);
                m_sparse_rmq.set_vector(&m_min_excess);
            }
            util::swap_support(m_sparse_rmq, rm.m_sparse_rmq,
                               &m_min_excess, &(rm.min_excess));
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
#ifdef MEASURE_TIMINGS
            if (t_block_size > 0) Stats::instance().addToTotal("Range",(r-l+1));
#endif
            TIME_MEASURE(size_type i = m_rank_select.select(l+2)-1; , "Select")

            num_queries++;
            if (r - l < 64) {
                value_type data = m_gct_bp.get_int(i+1);
                uint64_t one_cnt = bits::cnt(data)-1;
                if (l + one_cnt >= r) {
                    size_t cur_select = l+1, min_idx = i;
                    int cur_excess = ((l+1) << 1)-(i+1), min_excess = cur_excess;
                    for (size_t k = i+1; k <= i+65; ++k) {
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
                    num_avoided_selects++;
                    return (min_excess+min_idx)>>1;
                }
            }

            TIME_MEASURE(size_type j = m_rank_select.select(r+2); ,"Select")

            TIME_MEASURE(size_type sparse_i = (i+t_super_block_size-1)/t_super_block_size;
                         size_type sparse_j = j/t_super_block_size; ,"Other")
            bit_vector::difference_type min_rel_ex = 0;
            if (sparse_i >= sparse_j) {
                TIME_MEASURE(sparse_i = sparse_i-(i != 0); ,"Other")
                size_type rmq_e = i;
                if (sparse_i == sparse_j) {
                    TIME_MEASURE(rmq_e = get_min_excess_idx(sparse_i); ,"min_excess_idx")
                } else {
                    TIME_MEASURE(size_type rmq_e1 = get_min_excess_idx(sparse_i); ,"min_excess_idx")
                    TIME_MEASURE(size_type rmq_e2 = get_min_excess_idx(sparse_j); ,"min_excess_idx")
                    TIME_MEASURE(rmq_e = (m_min_excess[sparse_i] < m_min_excess[sparse_j] ? rmq_e1 : rmq_e2); , "min_excess")
                }
                if (rmq_e < i || rmq_e > j) {
                    TIME_MEASURE(rmq_e = near_rmq(m_gct_bp,i,j-1,min_rel_ex); ,"Scan");
//                   std::cout << i << " " << rmq_e << " " << j << " " << m_rank_select.rank(i+1)-1 << " " << m_rank_select.rank(rmq_e+1)-1 << " " << m_rank_select.rank(j+1)-2 << std::endl;
                }
                RETURN_TIME_MEASURE(size_type,m_rank_select.rank(rmq_e+1)-1; ,"Rank",t_block_size)
            } else {
                size_type rmq_sparse_idx = 0;
                if (t_block_size > 0) {
                    rmq_sparse_idx = (*m_rmq_recursive)(sparse_i,sparse_j-1);
                } else {
                    TIME_MEASURE(rmq_sparse_idx = m_sparse_rmq(sparse_i,sparse_j-1); ,"Sparse_RMQ")
                }
                TIME_MEASURE(size_type rmq_sparse = get_min_excess_idx(rmq_sparse_idx); ,"min_excess_idx")
                TIME_MEASURE(i_value_type rmq_sparse_ex = m_min_excess[rmq_sparse_idx]; ,"min_excess")
                TIME_MEASURE(size_type rmq_e1 = get_min_excess_idx(sparse_i-(i != 0)); ,"min_excess_idx")
                TIME_MEASURE(i_value_type rmq_e1_ex = m_min_excess[sparse_i-(i != 0)]; ,"min_excess")
                if (rmq_e1_ex < rmq_sparse_ex && rmq_e1 < i) {
                    TIME_MEASURE(rmq_e1 = near_rmq(m_gct_bp,i,t_super_block_size*sparse_i,min_rel_ex); ,"Scan")
                    TIME_MEASURE(rmq_e1_ex = m_rank_select.excess(rmq_e1); ,"Rank")
                }

                TIME_MEASURE(size_type rmq_e2 = get_min_excess_idx(sparse_j); ,"min_excess_idx");
                TIME_MEASURE(i_value_type rmq_e2_ex = m_min_excess[sparse_j]; ,"min_excess")
                if (rmq_e2_ex <= rmq_sparse_ex && rmq_e2 > j) {
                    TIME_MEASURE(rmq_e2 = near_rmq(m_gct_bp,t_super_block_size*sparse_j,j,min_rel_ex); ,"Scan")
                    TIME_MEASURE(rmq_e2_ex = m_rank_select.excess(rmq_e2); ,"Rank")
                }
                TIME_MEASURE(auto rmq_min = min_ex(rmq_e1,rmq_sparse,rmq_e2,
                                                   rmq_e1_ex,rmq_sparse_ex,rmq_e2_ex);
                             size_type rmq_e = rmq_min.first;
                             int_vector<>::value_type rmq_e_ex = rmq_min.second; ,"Other")
                RETURN_TIME_MEASURE(size_type,((rmq_e_ex+rmq_e)>>1); ,"Other",t_block_size)
            }
        }


        size_type size()const {
            return (m_gct_bp.size()-2)/2;
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_gct_bp.serialize(out, child, "gct_bp");
            written_bytes += m_rank_select.serialize(out, child, "rank_select_bp");
            written_bytes += m_min_excess.serialize(out, child, "min_excess");
            written_bytes += m_min_excess_idx.serialize(out, child, "min_excess_idx");
            if (t_block_size > 0) {
                written_bytes += m_rmq_recursive->serialize(out, child, "rmq_recursive");
            } else {
                written_bytes += m_sparse_rmq.serialize(out, child, "sparse_rmq");
            }
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in) {
            m_gct_bp.load(in);
            m_rank_select.load(in, &m_gct_bp);
            m_min_excess.load(in);
            m_min_excess_idx.load(in);
            if (t_block_size > 0) {
                m_rmq_recursive = new rmq_succinct_rec_old<t_block_size,0>();
                m_rmq_recursive->load(in);
            } else {
                m_sparse_rmq.load(in,&m_min_excess);
            }
        }
};

} // end namespace sdsl
#endif
