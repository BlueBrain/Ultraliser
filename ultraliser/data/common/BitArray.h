/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#ifndef ULTRALISER_DATA_BIT_ARRAY_H
#define ULTRALISER_DATA_BIT_ARRAY_H

#include <common/Common.h>

namespace Ultraliser
{

// Forward Declaration
class BitArray;

class BitArrayIndex
{
public:

    /**
     * @brief BitArrayIndex
     * This is the BitArrayIndex constructor.
     * It stores a pointer to the bit array and the bit index.
     * Pointer to bit array and bit index are stored.
     *
     * @param array
     * Pointer to bit array.
     *
     * @param index
     * Index of bit in array.
     */
    BitArrayIndex(BitArray *array, const uint64_t index);

public:
    /**
     * @brief operator =
     * Overload of the = operator.  Sets the bit array bit to the value of src.
     * Bit pointed to by this object is set to the value of source.
     * @param src Bit value.
     */
    void operator=(const bool src);

private:

    /**
     * @brief _bitArray
     * Array index applies to.
     */
    BitArray* _bitArray;

    /**
     * @brief _index
     * Index of bit in array.
     */
    uint64_t _index;
};

class BitArray
{
public:

    /**
     * @brief BitArray
     * This is the BitArray constructor. It reserves memory for the vector
     * storing the array. It Allocates vectory for array bits.
     * @param numBits Number of bits in the array.
     */
    BitArray(const uint64_t numBits);

    /**
     * @brief BitArray
     * This is the bit_array_c constructor.  It copies the for contents of a
     * vector of uint8_t into the BitArray.
     * @param array Vector to be copied.
     * @param numBits Number of bits in the array.
     */
    BitArray(uint8_t* array, const uint64_t numBits);

    ~BitArray() { delete[] _data; }

public:

    uint8_t* getData()
    {
        return _data;
    }

    /**
     * @brief print
     * Dumps the conents of a bit array to stdout.
     * The format of the dump is a series of bytes represented in hexadecimal.
     * @param outStream Stream to write to.
     */
    void print(std::ostream &outStream);

    /**
     * @brief writeBIN
     * @param prefix
     */
    void writeBIN(const std::string &prefix);

    /**
     * @brief size
     * @return
     */
    uint64_t size() const { return _numBits; }

    /**
     * @brief getNumberBytes
     * Total number of bytes in the array.
     * @return
     */
    uint64_t getNumberBytes() const { return _numBytes; }

    /**
     * @brief setAll
     * Sets every bit in the bit array to 1.
     * This method uses uint8_t_MAX to determine what it means to set all bits
     * in an uint8_t, so it is crucial that the machine implementation of uint8_t
     * utilizes all of the bits in the memory allocated for an uint8_t.
     * Each of the bits used in the bit array are set to 1.
     * Unused (spare) bits are set to 0.
     */
    void setAll();

    /**
     * @brief clearAll
     * Sets every bit in the bit array to 0.
     * Each of the bits in the bit array are set to 0.
     */
    void clearAll();

    /**
     * @brief setBit
     * Sets a bit in the bit array to 1.
     * The specified bit will be set to 1.
     * @param bit The number of the bit to set.
     */
    void setBit(const uint64_t bit);

    /**
     * @brief clearBit
     * Sets a bit in the bit array to 0.
     * The specified bit will be set to 0.
     * @param bit The number of the bit to clear.
     */
    void clearBit(const uint64_t bit);

    /**
     * @brief bit
     * @param index
     * @return
     */
    bool bit(const uint64_t index) const;

    /**
     * @brief getByte
     * @param index
     * @return
     */
    uint8_t getByte(const uint64_t index) const;

    /**
     * @brief addByte
     * @param index
     * @param byte
     */
    void addByte(const uint64_t index, const uint8_t byte);

public:

    /**
     * @brief operator ()
     * Overload of the () operator.
     * This method approximates array indices used for assignment.
     * It returns a BitArrayIndex which includes an = method used to
     * set bit values.
     * @param bit Index of array bit.
     * @return
     */
    BitArrayIndex operator()(const uint64_t bit);

    /**
     * @brief operator []
     * Overload of the [] operator.
     * This method returns the value of a bit in the bit array.
     * @param bit Index of array bit.
     * @return The value of the specified bit.
     */
    bool operator[](const uint64_t bit) const;

    /**
     * @brief operator ==
     * Overload of the == operator.
     * @param other BitArray to compare.
     * @return True if this == other.  Otherwise false.
     */
    bool operator==(const BitArray &other) const;

    /**
     * @brief operator !=
     * Overload of the != operator.
     * @param other BitArray to compare.
     * @return True if this != other.  Otherwise false.
     */
    bool operator!=(const BitArray &other) const;

    /**
     * @brief operator <
     * Overload of the < operator.
     * @param other BitArray to compare.
     * @return True if this < other.  Otherwise
     */
    bool operator<(const BitArray &other) const;

    /**
     * @brief operator <=
     * Overload of the <= operator.
     * @param other BitArray to compare.
     * @return True if this <= other.  Otherwise false.
     */
    bool operator<=(const BitArray &other) const;

    /**
     * @brief operator >
     * Overload of the > operator.
     * @param other BitArray to compare.
     * @return True if this > other.  Otherwise false.
     */
    bool operator>(const BitArray &other) const;

    /**
     * @brief operator >=
     * Overload of the >= operator.
     * @param other BitArray to compare.
     * @return True if this >= other.  Otherwise false.
     */
    bool operator>=(const BitArray &other) const;

    /**
     * @brief operator &
     * Overload of the & operator. Performs a bitwise and between the source
     * array and this bit array.
     * @param other BitArray on righthand side of &.
     * @return Value of bitwise and of this and other.
     */
    BitArray operator&(const BitArray &other) const;

    /**
     * @brief operator ^
     * Overload of the ^ operator.  Performs a bitwise xor between the source
     * array and this bit array.
     * @param other BitArray on righthand side of ^.
     * @return Value of bitwise xor of this and other.
     */
    BitArray operator^(const BitArray &other) const;

    /**
     * @brief operator |
     * Overload of the | operator.  Performs a bitwise or between the source
     * array and this bit array.
     * @param other BitArray on righthand side of |.
     * @return Value of bitwise or of this and other.
     */
    BitArray operator|(const BitArray &other) const;

    /**
     * @brief operator ~
     * Overload of the ~ operator.  Negates all non-spare bits in BitArray.
     * @return Value of this after bitwise not.
     */
    BitArray operator~() const;

    /**
     * @brief operator <<
     * Overload of the << operator. Performs a bitwise left shift of this
     * BitArray.
     * @param count The number of bits to shift left.
     * @return Result of bitwise left shift.
     */

    BitArray operator<<(const uint64_t count) const;

    /**
     * @brief operator >>
     * Overload of the >> operator.  Performs a bitwise right shift of this
     * BitArray.
     * @param count The number of bits to shift right.
     * @return Result of bitwise right shift.
     */
    BitArray operator>>(const uint64_t count) const;

    /**
     * @brief operator ++
     * Overload of the ++ operator. Increments the contents of a BitArray.
     * Overflows cause rollover.
     * Bit array contents are incremented.
     * @return Reference to this array after increment.
     */
    BitArray& operator++();

    /**
     * @brief operator ++
     * Overload of the ++ operator. Increments the contents of a BitArray.
     * Overflows cause rollover.
     * Bit array contents are incremented.
     * @param dummy Needed for postfix increment.
     * @return Reference to this array after increment.
     */
    BitArray& operator++(int dummy);

    /**
     * @brief operator --
     * Overload of the -- operator.
     * Decrements the contents of a BitArray.
     * Underflows cause rollover.
     * Bit array contents are decremented.
     * @return
     */
    BitArray& operator--();

    /**
     * @brief operator --
     * Overload of the -- operator.
     * Decrements the contents of a BitArray.
     * Underflows cause rollover.
     * Bit array contents are decremented.
     * @param dummy Needed for postfix decrement.
     * @return
     */
    BitArray& operator--(int dummy);

    /**
     * @brief operator =
     * Overload of the = operator.
     * Copies source contents into this BitArray.
     * Source BitArray contents are copied into this array.
     * @param src Source BitArray.
     * @return Reference to this array after copy.
     */
    BitArray& operator=(const BitArray &src);

    /**
     * @brief operator &=
     * Overload of the &= operator. Performs a bitwise and between the source
     * array and this bit array.  This BitArray will contain the result.
     * Results of bitwise and are stored in this array
     * @param src Source BitArray.
     * @return Reference to this array after and.
     */
    BitArray& operator&=(const BitArray &src);

    /**
     * @brief operator ^=
     * Overload of the ^= operator. Performs a bitwise xor between the source
     * array and this BitArray.
     * This bit array will contain the result.
     * Results of bitwise xor are stored in this array.
     * @param src Source BitArray.
     * @return Reference to this array after xor.
     */
    BitArray& operator^=(const BitArray &src);

    /**
     * @brief operator |=
     * Overload of the |= operator.
     * Performs a bitwise or between the source array and this BitArray.
     * This BitArray will contain the result.
     * Results of bitwise or are stored in this array.
     * @param src Source BitArray.
     * @return Reference to this array after or
     */
    BitArray& operator|=(const BitArray &src);

    /**
     * @brief Not
     * Negates the array, negate (~=).
     * Negates all non-spare bits in bit array.
     * Contents of bit array are negated.  Any spare bits are left at 0.
     * @return Reference to this array after shift.
     */
    BitArray& Not();

    /**
     * @brief operator <<=
     * Overload of the <<= operator.  Performs a left shift on this bit array.
     * This bit array will contain the result.
     * Results of the shifts are stored in this array.
     * @param shifts Number of bit positions to shift.
     * @return Reference to this array after shift.
     */
    BitArray& operator<<=(const uint64_t shifts);

    /**
     * @brief operator >>
     * Overload of the >>= operator.  Performs a right shift on this bit array.
     * This bit array will contain the result.
     * Results of the shifts are stored in this array.
     * @param shifts Number of bit positions to shift.
     * @return Reference to this array after shift.
     */
    BitArray& operator>>=(const uint64_t shifts);

protected:

    /**
     * @brief _numBits
     * Number of bits in the array.
     */
    uint64_t _numBits;

    /**
     * @brief _numBytes
     */
    uint64_t _numBytes;

    /**
     * @brief data_
     * Vector of characters.
     */
    uint8_t *_data;
};

}
#endif // ULTRALISER_DATA_BIT_ARRAY_H
