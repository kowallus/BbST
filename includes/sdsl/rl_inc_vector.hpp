/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

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
/*! \file rl_inc_vector.hpp
   \brief rl_inc_vector.hpp contains the sdsl::rl_inc_vector class.
   \author Simon Gog
*/
#ifndef SDSL_RL_INC_VECTOR
#define SDSL_RL_INC_VECTOR

#include "int_vector.hpp"
#include "dac_vector.hpp"
#include "bit_vectors.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "coder.hpp"
#include "iterators.hpp"

#include <set>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A generic immutable space-saving vector class for unsigned integers.
/*! A vector v is stored more space-efficiently by self-delimiting coding
 *  the deltas v[i+1]-v[i] (v[-1]:=0). Space of the structure and random
 *  access time to it can be controlled by a sampling parameter t_dens.
 *
 *  \tparam t_coder  Self-delimiting coder.
 *  \tparam t_dens   Every t_dens-th element of v is sampled.
 *  \tparam t_width  Width of the int_vector used to store the samples and pointers.
 *  This class is a parameter of csa_sada.
 * @ingroup int_vector
 */
template <class t_int_vector = dac_vector<>, class t_bit_vector = sd_vector<>,
          uint32_t t_dens = 8, uint8_t t_width = 0>
class rl_inc_vector
{
  private:
    //static_assert(t_dens > 1 , "rl_inc_vector: sample density must be larger than `1`");
  public:
    typedef uint64_t value_type;
    typedef random_access_const_iterator<rl_inc_vector> iterator;
    typedef iterator const_iterator;
    typedef ptrdiff_t difference_type;
    typedef int_vector<>::size_type size_type;
    typedef iv_tag index_category;
    typedef typename t_bit_vector::rank_1_type rank_support;
    typedef typename t_bit_vector::select_1_type select_support;

    static const uint32_t sample_dens = t_dens;

  private:
    size_type m_size = 0; // number of vector elements

    t_int_vector m_differences;

    t_bit_vector m_samples;
    select_support m_samples_select;

    t_bit_vector m_run_marker;
    rank_support m_run_rank;
    select_support m_run_select;

    size_type m_run_count;

    void clear()
    {
        m_size = 0;
        m_run_count = 0;
    }

    inline size_type run_length(size_type run_idx) const
    {
        return m_run_select(run_idx + 2) - m_run_select(run_idx + 1) - 1;
    }

    inline value_type encode_difference(const size_type run_idx, const value_type value_before) const
    {
        return m_differences[run_idx] + value_before;
    }

  public:
    rl_inc_vector() = default;
    rl_inc_vector(const rl_inc_vector &) = default;
    rl_inc_vector(rl_inc_vector &&) = default;
    rl_inc_vector &operator=(const rl_inc_vector &) = default;
    rl_inc_vector &operator=(rl_inc_vector &&) = default;

    //! Constructor for a Container of unsigned integers.
    /*! \param c A container of unsigned integers.
          */
    template <class Container>
    rl_inc_vector(const Container &c);

    //! Constructor for an int_vector_buffer of unsigned integers.
    /*
            \param v_buf A int_vector_buf.
        */
    template <uint8_t int_width>
    rl_inc_vector(int_vector_buffer<int_width> &v_buf);

    //! Default Destructor
    ~rl_inc_vector() {}

    //! The number of elements in the rl_inc_vector.
    size_type size() const
    {
        return m_size;
    }

    //! Return the largest size that this container can ever have.
    static size_type max_size()
    {
        return int_vector<>::max_size() / 2;
    }

    //!    Returns if the rl_inc_vector is empty.
    bool empty() const
    {
        return 0 == m_size;
    }

    //! Swap method for rl_inc_vector
    void swap(rl_inc_vector &v);

    //! Iterator that points to the first element of the rl_inc_vector.
    const const_iterator begin() const
    {
        return const_iterator(this, 0);
    }

    //! Iterator that points to the position after the last element of the rl_inc_vector.
    const const_iterator end() const
    {
        return const_iterator(this, this->m_size);
    }

    //! operator[]
    /*! \param i Index. \f$ i \in [0..size()-1]\f$.
         */
    value_type operator[](size_type i) const;

    //! Serialize the rl_inc_vector to a stream.
    /*! \param out Out stream to write the data structure.
            \return The number of written bytes.
         */
    size_type serialize(std::ostream &out, structure_tree_node *v = nullptr, std::string name = "") const;

    //! Load the rl_inc_vector from a stream.
    void load(std::istream &in);

    //! Returns the i-th sample of rl_inc_vector
    /*! \param i The index of the sample. 0 <= i < size()/get_sample_dens()
         *  \return The value of the i-th sample.
         */
    //value_type sample(const size_type i) const;

    uint32_t get_sample_dens() const
    {
        return t_dens;
    }

    /*!
         * \param i The index of the sample for which all values till the next sample should be decoded. 0 <= i < size()/get_sample_dens()
         * \param it A pointer to a uint64_t vector, whereto the values should be written
         */
    /*void get_inter_sampled_values(const size_type i, uint64_t* it)const
        {
            *(it++) = 0;
            if (i*t_dens + t_dens - 1 < size()) {
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], t_dens - 1, it);
            } else {
                assert(i*t_dens < size());
                t_coder::template decode<true, true>(m_z.data(), m_sample_vals_and_pointer[(i<<1)+1], size()-i*t_dens - 1, it);
            }
        };*/
};

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
inline typename rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::value_type
    rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::operator[](const size_type i) const
{
    size_type run_idx = m_run_rank(i + 1) - 1;
    size_type run_idx_select = m_run_select(run_idx + 1);
    size_type sample_idx = run_idx / t_dens;
    size_type cur_run_idx = sample_idx * t_dens;
    value_type val = m_samples_select(sample_idx+1);
    //std::cout << sample_idx << " " << val << std::endl;
    for (; cur_run_idx <= run_idx; ++cur_run_idx)
    {
        val = encode_difference(cur_run_idx, val);
        if (cur_run_idx < run_idx)
        {
            val += run_length(cur_run_idx);
        }
    }
    val += i - run_idx_select;
    return val;
}

/*template<class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
inline typename rl_inc_vector<t_coder, t_dens,t_width>::value_type rl_inc_vector<t_coder, t_dens,t_width>::sample(const size_type i)const
{
    assert(i*get_sample_dens()+1 != 0);
    assert(i*get_sample_dens() < m_size);
    return m_sample_vals_and_pointer[i<<1];
}*/

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
void rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::swap(rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width> &v)
{
    if (this != &v)
    {
        std::swap(m_size, v.m_size);
        std::swap(m_run_count, v.m_run_count);
        m_differences.swap(v.m_differences);
        m_samples.swap(v.m_samples);
        util::swap_support(m_samples_select, v.m_samples_select,
                           &m_samples, &(v.m_samples));
        m_run_marker.swap(v.m_run_marker);
        util::swap_support(m_run_rank, v.m_run_rank,
                           &m_run_marker, &(v.m_run_marker));
        util::swap_support(m_run_select, v.m_run_select,
                           &m_run_marker, &(v.m_run_marker));
    }
}

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
template <class Container>
rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::rl_inc_vector(const Container &c)
{

    // clear bit_vectors
    clear();

    m_size = c.size();

    int_vector<> differences = int_vector<>(m_size, 0);
    bit_vector samples = bit_vector(c[m_size-1]+1, 0);
    samples[0] = 1;

    bit_vector run_marker = bit_vector(m_size + 1, 0);
    run_marker[0] = 1;
    run_marker[m_size] = 1;

    m_run_count = 0;
    differences[m_run_count++] = c[0];
    for (size_type i = 1; i < m_size; ++i)
    {
        if (c[i] != c[i - 1] + 1)
        {
            run_marker[i] = 1;
            if (c[i] >= c[i - 1])
            {
                differences[m_run_count] = c[i] - c[i - 1];
            }
            else
            {
                differences[m_run_count] = std::numeric_limits<value_type>::max() - (c[i - 1] - c[i]) + 1;
            }
            if (m_run_count % t_dens == 0)
            {
                samples[c[i - 1]] = 1;
            }
            m_run_count++;
        }
    }
    differences.resize(m_run_count);

    m_differences = t_int_vector(differences);
    m_run_marker = t_bit_vector(run_marker);
    m_samples = t_bit_vector(samples);

    m_run_rank = rank_support(&m_run_marker);
    m_run_select = select_support(&m_run_marker);
    m_samples_select = select_support(&m_samples);
}

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
template <uint8_t int_width>
rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::rl_inc_vector(int_vector_buffer<int_width> &v_buf)
{
    // clear bit_vectors
    clear();

    m_size = v_buf.size();

    int_vector<> differences = int_vector<>(m_size, 0);
    bit_vector samples = bit_vector(v_buf[m_size - 1] + 1, 0);
    samples[0] = 1;

    bit_vector run_marker = bit_vector(m_size + 1, 0);
    run_marker[0] = 1;
    run_marker[m_size] = 1;

    m_run_count = 0;
    differences[m_run_count++] = v_buf[0];
    for (size_type i = 1; i < m_size; ++i)
    {
        if (v_buf[i] != v_buf[i - 1] + 1)
        {
            run_marker[i] = 1;
            if (v_buf[i] >= v_buf[i - 1])
            {
                differences[m_run_count] = v_buf[i] - v_buf[i - 1];
            }
            else
            {
                differences[m_run_count] = std::numeric_limits<value_type>::max() - (v_buf[i - 1] - v_buf[i]) + 1;
            }
            if (m_run_count % t_dens == 0)
            {
                samples[v_buf[i - 1]] = 1;
            }
            m_run_count++;
        }
    }

    differences.resize(m_run_count);

    m_differences = t_int_vector(differences);
    m_run_marker = t_bit_vector(run_marker);
    m_samples = t_bit_vector(samples);

    m_run_rank = rank_support(&m_run_marker);
    m_run_select = select_support(&m_run_marker);
    m_samples_select = select_support(&m_samples);
}

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
rl_inc_vector<>::size_type rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::serialize(std::ostream &out, structure_tree_node *v, std::string name) const
{
    structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "size");
    written_bytes += write_member(m_run_count, out, child, "number runs");
    written_bytes += m_differences.serialize(out, child, "differences");
    written_bytes += m_samples.serialize(out, child, "samples");
    written_bytes += m_samples_select.serialize(out, child, "samples select");
    written_bytes += m_run_marker.serialize(out, child, "run marker");
    written_bytes += m_run_rank.serialize(out, child, "run marker rank");
    written_bytes += m_run_select.serialize(out, child, "run marker select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

template <class t_int_vector, class t_bit_vector, uint32_t t_dens, uint8_t t_width>
void rl_inc_vector<t_int_vector, t_bit_vector, t_dens, t_width>::load(std::istream &in)
{
    read_member(m_size, in);
    read_member(m_run_count, in);
    m_differences.load(in);
    m_samples.load(in);
    m_samples_select.load(in, &m_samples);
    m_run_marker.load(in);
    m_run_rank.load(in, &m_run_marker);
    m_run_select.load(in, &m_run_marker);
}

} // end namespace sdsl
#endif
