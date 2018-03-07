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
/*! \file rl_bitenc_vector.hpp
   \brief rl_bitenc_vector.hpp contains the sdsl::rl_bitenc_vector class.
   \author Simon Gog
*/
#ifndef SDSL_rl_bitenc_vector
#define SDSL_rl_bitenc_vector

#include "int_vector.hpp"
#include "bit_vectors.hpp"
#include "coder.hpp"
#include "iterators.hpp"


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
class rl_bitenc_vector
{
    public:
        typedef uint64_t                                 value_type;
        typedef random_access_const_iterator<rl_bitenc_vector>  iterator;
        typedef iterator                                 const_iterator;
        typedef const value_type                         reference;
        typedef const value_type                         const_reference;
        typedef const value_type*                        const_pointer;
        typedef ptrdiff_t                                difference_type;
        typedef int_vector<>::size_type                  size_type;
        typedef iv_tag                                   index_category;
        typedef rl_bitenc_vector                         rl_vec_type;
        typedef typename rl_bit_vector<>::rank_1_type    rank_support;
        typedef typename rl_bit_vector<>::select_1_type  select_support;

    private:
        size_type         m_size = 0;                // number of vector elements

        rl_bit_vector<> m_rl_bit;
        select_support m_rl_bit_select;

        void clear()
        {
            m_size = 0;
        }

    public:
        rl_bitenc_vector() = default;
        rl_bitenc_vector(const rl_bitenc_vector&) = default;
        rl_bitenc_vector(rl_bitenc_vector&&) = default;
        rl_bitenc_vector& operator=(const rl_bitenc_vector&) = default;
        rl_bitenc_vector& operator=(rl_bitenc_vector&&) = default;

        //! Constructor for a Container of unsigned integers.
        /*! \param c A container of unsigned integers.
          */
        template<class Container>
        rl_bitenc_vector(const Container& c);

        //! Constructor for an int_vector_buffer of unsigned integers.
        /*
            \param v_buf A int_vector_buf.
        */
        template<uint8_t int_width>
        rl_bitenc_vector(int_vector_buffer<int_width>& v_buf);

        //! Default Destructor
        ~rl_bitenc_vector() { }

        //! The number of elements in the rl_bitenc_vector.
        size_type size()const
        {
            return m_size;
        }

        //! Return the largest size that this container can ever have.
        static size_type max_size()
        {
            return int_vector<>::max_size()/2;
        }

        //!    Returns if the rl_bitenc_vector is empty.
        bool empty() const
        {
            return 0==m_size;
        }

        //! Swap method for rl_bitenc_vector
        void swap(rl_bitenc_vector& v);

        //! Iterator that points to the first element of the rl_bitenc_vector.
        const const_iterator begin()const
        {
            return const_iterator(this, 0);
        }

        //! Iterator that points to the position after the last element of the rl_bitenc_vector.
        const const_iterator end()const
        {
            return const_iterator(this, this->m_size);
        }

        //! operator[]
        /*! \param i Index. \f$ i \in [0..size()-1]\f$.
         */
        value_type operator[](size_type i)const;

        //! Serialize the rl_bitenc_vector to a stream.
        /*! \param out Out stream to write the data structure.
            \return The number of written bytes.
         */
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const;

        //! Load the rl_bitenc_vector from a stream.
        void load(std::istream& in);

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

inline typename rl_bitenc_vector::value_type rl_bitenc_vector::operator[](const size_type i)const
{
    assert(i+1 != 0);
    assert(i < m_size);
    
    value_type value = m_rl_bit_select(i+1);
    return value;
}

void rl_bitenc_vector::swap(rl_bitenc_vector& v)
{
    if (this != &v) { // if v and _this_ are not the same object
        std::swap(m_size, v.m_size);
        std::swap(m_rl_bit, v.m_rl_bit);
        util::swap_support(m_rl_bit_select, v.m_rl_bit_select,
                    &m_rl_bit, &(v.m_rl_bit));
    }
}

template<class Container>
rl_bitenc_vector::rl_bitenc_vector(const Container& c)
{
    // clear bit_vectors
    clear();

    if (c.empty())  // if c is empty there is nothing to do...
        return;

    m_rl_bit = rl_bit_vector<>(c.begin(),c.end());
    util::init_support(m_rl_bit_select, &m_rl_bit);
}

template<uint8_t int_width>
rl_bitenc_vector::rl_bitenc_vector(int_vector_buffer<int_width>& v_buf)
{

}

rl_bitenc_vector::size_type rl_bitenc_vector::serialize(std::ostream& out, structure_tree_node* v, std::string name)const
{
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "size");
    written_bytes += m_rl_bit.serialize(out, child, "rl bit vector");
    written_bytes += m_rl_bit_select.serialize(out, child, "rl bit vector select");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void rl_bitenc_vector::load(std::istream& in)
{
    read_member(m_size, in);
    m_rl_bit.load(in);
    m_rl_bit_select.load(in, &m_rl_bit);
}

} // end namespace sdsl
#endif
