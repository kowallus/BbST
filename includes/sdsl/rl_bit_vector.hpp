/* sdsl - succinct data structures library
    Copyright (C) 2012-2014 Simon Gog
    Copyright (C) 2015 Genome Research Ltd.

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
/*!\file rl_bit_vector.hpp
   \brief rl_bit_vector.hpp contains the sdsl::rl_bit_vector class, and
          classes which support rank and select for rl_bit_vector.
   \author Simon Gog, Jouni Siren
*/
#ifndef INCLUDED_SDSL_RL_BIT_VECTOR
#define INCLUDED_SDSL_RL_BIT_VECTOR

#include "util.hpp"
#include "iterators.hpp"

//! Namespace for the succinct data structure library
namespace sdsl
{

template<uint8_t t_b          = 1,
         class t_bit_vector   = bit_vector>
class rank_support_rl; 

template<uint8_t t_b          = 1,
         class t_bit_vector   = bit_vector>
class select_support_rl;

//!
/*!
 * 
 */
template<class t_bit_vector = sd_vector<>>
class rl_bit_vector
{
    public:
        typedef bit_vector::size_type                       size_type;
        typedef size_type                                   value_type;
        typedef bit_vector::difference_type                 difference_type;
        typedef random_access_const_iterator<rl_bit_vector> iterator;
        typedef iterator                                    const_iterator;
        typedef t_bit_vector                                bit_vector_type;
        typedef typename t_bit_vector::rank_1_type          rank_support;
        typedef typename t_bit_vector::select_1_type        select_support;

        typedef rank_support_rl<0, t_bit_vector> rank_0_type;
        typedef rank_support_rl<1, t_bit_vector> rank_1_type;
        typedef select_support_rl<0, t_bit_vector> select_0_type;
        typedef select_support_rl<1, t_bit_vector> select_1_type;
    private:
    
        size_type m_size = 0;  // length of the original bit vector

        bit_vector_type m_b1;
        bit_vector_type m_b_rl; 
        rank_support    m_b1_rank_1;
        rank_support    m_b_rl_rank_1;
        select_support  m_b1_select_1;
        select_support  m_b_rl_select_1;

        void copy(const rl_bit_vector& v)
        {
            m_size = v.m_size;
            m_b1 = v.m_b1;
            m_b_rl = v.m_b_rl;
            m_b1_rank_1 = v.m_b1_rank_1;
            m_b1_select_1 = v.m_b1_select_1;
            m_b_rl_rank_1 = v.m_b_rl_rank_1;
            m_b_rl_select_1 = v.m_b_rl_select_1;
        }

    public:
        const bit_vector_type& b1            = m_b1;
        const bit_vector_type& b_rl          = m_b_rl; 
        const rank_support&    b1_rank_1     = m_b1_rank_1;
        const rank_support&    b_rl_rank_1   = m_b_rl_rank_1;
        const select_support&  b1_select_1   = m_b1_select_1;
        const select_support&  b_rl_select_1 = m_b_rl_select_1;

        rl_bit_vector() { }

        rl_bit_vector(const rl_bit_vector& rl)
        {
            copy(rl);
        }

        rl_bit_vector(rl_bit_vector&& rl)
        {
            *this = std::move(rl);
        }

        rl_bit_vector(const bit_vector& bv)
        {
            if(bv.size() > 0) {
                m_size = bv.size();
                size_type total_runs_length = bv[0];

                //(1) Setup bit vector with run length heads
                bit_vector tmp_b1(m_size,0);
                tmp_b1[0] = bv[0];
                for(size_type i = 1; i < m_size; ++i) {
                    if(bv[i] == 1 && bv[i-1] == 0) {
                        tmp_b1[i] = 1; 
                    }
                    total_runs_length += bv[i];
                }
                util::assign(m_b1, tmp_b1);
                util::init_support(m_b1_rank_1, &m_b1);
                util::init_support(m_b1_select_1, &m_b1);

                //(2) Setup bit vector with run lengths in unary encoding
                bit_vector tmp_b_rl(total_runs_length,0);
                size_t idx = 0;
                for(size_type i = 0; i < m_size-1; ++i) {
                    if(bv[i] == 1 && bv[i+1] == 0) {
                        tmp_b_rl[idx] = 1;
                    }
                    idx += bv[i];
                }
                tmp_b_rl[idx] = bv[m_size-1];
                util::assign(m_b_rl, tmp_b_rl);
                util::init_support(m_b_rl_rank_1, &m_b_rl);
                util::init_support(m_b_rl_select_1, &m_b_rl);
            }
        }

        //TODO(heuer): Write explicit construction from iterator
        template<class t_itr>
        rl_bit_vector(const t_itr begin,const t_itr end)
        {
            if (begin == end) {
                return;
            }
            if (! is_sorted(begin,end)) {
                throw std::runtime_error("rl_bit_vector: source list is not sorted.");
            }

            size_type m = std::distance(begin,end);
            m_size = *(end-1)+1;
            int_vector<> run_heads(m);
            int_vector<> run_length(m);

            t_itr cur = begin;
            run_heads[0] = *cur;
            value_type last_value = *(cur++);
            size_type run_start = 0, num_runs = 1;
            for(size_type i = 1; cur != end; ++cur, ++i) {
                value_type value = *cur;
                if(value != last_value + 1) {
                    run_length[num_runs-1] = i-1;
                    run_heads[num_runs++] = value;
                }
                last_value = value;
            }
            run_length[num_runs-1] = m-1;
            run_heads.resize(num_runs);
            run_length.resize(num_runs);

            m_b1 = t_bit_vector(run_heads.begin(),run_heads.end());
            m_b_rl = t_bit_vector(run_length.begin(),run_length.end());

            util::init_support(m_b1_rank_1, &m_b1);
            util::init_support(m_b1_select_1, &m_b1);
            util::init_support(m_b_rl_rank_1, &m_b_rl);
            util::init_support(m_b_rl_select_1, &m_b_rl);
        }


        //! 
        /*! 
        */
        value_type operator[](size_type i)const
        {
            size_type run_idx = m_b1_rank_1(i+1);
            if(run_idx == 0) return 0;
            size_type run_start_pos = m_b1_select_1(run_idx);
            size_type run_length = get_run_length(run_idx);
            return run_start_pos + run_length > i;
        }

        //! Get the integer value of the binary string of length len starting at position idx.
        /*! \param idx Starting index of the binary representation of the integer.
         *  \param len Length of the binary representation of the integer. Default value is 64.
         *  \returns The integer value of the binary string of length len starting at position idx.
         *
         *  \pre idx+len-1 in [0..size()-1]
         *  \pre len in [1..64]
         */
        uint64_t get_int(size_type idx, const uint8_t len=64) const
        {
            uint64_t res = 0;
            size_type cur_idx = 0;
            size_type cur_run = m_b1_rank_1(idx+1);
            size_type len_64_t = std::min(static_cast<size_type>(len), m_size - idx + 1);

            if((*this)[idx]) {
                size_type run_start = m_b1_select_1(cur_run);
                size_type run_length = get_run_length(cur_run);
                size_type partial_run_length = std::min(run_start + run_length - idx, len_64_t);
                res |= ((1ULL << partial_run_length) - 1);
                cur_idx += partial_run_length;
            } 
                
            size_type run_start = m_b1_select_1(++cur_run);
            cur_idx += (run_start - (idx + cur_idx));

            while(cur_idx < len_64_t) {
                size_type run_length = std::min(get_run_length(cur_run), len_64_t - cur_idx);
                res |= (((1ULL << run_length) - 1) << cur_idx);
                size_type next_run_start_pos = m_b1_select_1(++cur_run);
                cur_idx += (next_run_start_pos - (idx + cur_idx));
            }
            return res;
        }

        //! Swap method
        void swap(rl_bit_vector& v)
        {
            if (this != &v) {
                std::swap(m_size, v.m_size);
                m_b1.swap(v.m_b1);
                m_b_rl.swap(v.m_b_rl);
                util::swap_support(m_b1_rank_1, v.m_b1_rank_1, &m_b1, &v.m_b1);
                util::swap_support(m_b1_select_1, v.m_b1_select_1, &m_b1, &v.m_b1);
                util::swap_support(m_b_rl_rank_1, v.m_b_rl_rank_1, &m_b_rl, &v.m_b_rl);
                util::swap_support(m_b_rl_select_1, v.m_b_rl_select_1, &m_b_rl, &v.m_b_rl);
            }
        }

        //! Returns the size of the original bit vector.
        size_type size()const
        {
            return m_size;
        }

        rl_bit_vector& operator=(const rl_bit_vector& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        rl_bit_vector& operator=(rl_bit_vector&& v)
        {
            if (this != &v) {
                m_size = v.m_size;
                m_b1 = std::move(v.m_b1);
                m_b_rl = std::move(v.m_b_rl);
                m_b1_rank_1 = std::move(v.m_b1_rank_1);
                m_b1_rank_1.set_vector(&m_b1);
                m_b1_select_1 = std::move(v.m_b1_select_1);
                m_b1_select_1.set_vector(&m_b1);
                m_b_rl_rank_1 = std::move(v.m_b_rl_rank_1);
                m_b_rl_rank_1.set_vector(&m_b_rl);
                m_b_rl_select_1 = std::move(v.m_b_rl_select_1);
                m_b_rl_select_1.set_vector(&m_b_rl);
            }
            return *this;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "size");
            written_bytes += m_b1.serialize(out, child, "b1");
            written_bytes += m_b_rl.serialize(out, child, "b_rl");
            written_bytes += m_b1_rank_1.serialize(out, child, "m_b1_rank");
            written_bytes += m_b1_select_1.serialize(out, child, "m_b1_select");
            written_bytes += m_b_rl_rank_1.serialize(out, child, "m_b_rl_rank");
            written_bytes += m_b_rl_select_1.serialize(out, child, "m_b_rl_select");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in)
        {
            read_member(m_size,in);
            m_b1.load(in);
            m_b_rl.load(in);
            m_b1_rank_1.load(in, &m_b1);
            m_b1_select_1.load(in, &m_b1);
            m_b_rl_rank_1.load(in, &m_b_rl);
            m_b_rl_select_1.load(in, &m_b_rl);
        }

        inline size_type get_run_length(const size_type run_idx) const {
            size_type i = run_idx == 1 ? 0 : m_b_rl_select_1(run_idx - 1) + 1;
            size_type j = m_b_rl_select_1(run_idx);
            return j - i + 1;
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }
};


template<uint8_t t_b>
struct rank_support_rl_trait {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type)
    {
        return r;
    }
};

template<>
struct rank_support_rl_trait<0> {
    typedef bit_vector::size_type size_type;
    static size_type adjust_rank(size_type r, size_type n)
    {
        return n - r;
    }
};


//! Rank data structure for rl_bit_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_bit_vector    Type of the bitvector used for the run head and
 *                          run length bitvector.
 */
template<uint8_t t_b, class t_bit_vector>
class rank_support_rl
{
        static_assert(t_b == 1u or t_b == 0u , "rank_support_rl: bit pattern must be `0` or `1`");
    public:
        typedef bit_vector::size_type size_type;
        typedef rl_bit_vector<t_bit_vector> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;

    public:

        explicit rank_support_rl(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type rank(size_type i)const
        {
            assert(m_v != nullptr);
            assert(i <= m_v->size());
            if(m_v->size() == 0) return 0;
            size_type r = m_v->b1_rank_1(i);
            if(r == 0) return t_b ? 0 : i;
            size_type j = m_v->b1_select_1(r);
            size_type rank_j = (r == 1 ? 0 : m_v->b_rl_select_1(r-1) + 1);
            size_type k = m_v->b_rl_select_1(r) - rank_j + 1;
            size_type res = rank_j + (i-j >= k ? k : (i - j));
            return rank_support_rl_trait<t_b>::adjust_rank(res, i);
        }

        size_type operator()(size_type i)const
        {
            return rank(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        rank_support_rl& operator=(const rank_support_rl& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        void swap(rank_support_rl&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }
};


//! Select data structure for sd_vector
/*! \tparam t_b             Bit pattern.
 *  \tparam t_hi_bit_vector Type of the bitvector used for the unary decoded differences of
 *                          the high part of the positions of the 1s.
 *  \tparam t_select_1      Type of the select structure which is used to select ones in HI.
 *  \tparam t_select_0      Type of the select structure which is used to select zeros in HI.
 */
template<uint8_t t_b, class t_bit_vector>
class select_support_rl
{
    public:
        typedef bit_vector::size_type size_type;
        typedef rl_bit_vector<t_bit_vector> bit_vector_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;
    public:

        explicit select_support_rl(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const
        {
            size_type r = m_v->b_rl_rank_1(i-1) + 1;
            size_type k = (i-1) - (r == 1 ? 0 : m_v->b_rl_select_1(r-1) + 1);
            return m_v->b1_select_1(r) + k ;
        }

        size_type operator()(size_type i)const
        {
            return select(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        select_support_rl& operator=(const select_support_rl& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_support_rl&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }
};

template<class t_bit_vector>
class select_support_rl<0, t_bit_vector>
{
    public:
        typedef bit_vector::size_type size_type;
        typedef rl_bit_vector<t_bit_vector> bit_vector_type;
        typedef typename bit_vector_type::rank_0_type rank_support;
        enum { bit_pat = 0 };
        enum { bit_pat_len = (uint8_t)1 };
    private:

        const bit_vector_type* m_v;
        rank_support rank0;
    public:

        explicit select_support_rl(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const
        {
            size_type l = 0, r = m_v->size();
            while(l < r) {
                size_type m = (l+r)/2;
                size_type select0 = rank0(m);
                if(select0 < i) l = m+1;
                else r = m; 
            }
            return l-1;
        }

        size_type operator()(size_type i)const
        {
            return select(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
            rank0.set_vector(v);
        }

        select_support_rl& operator=(const select_support_rl& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_support_rl&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }
};


} // end namespace
#endif
