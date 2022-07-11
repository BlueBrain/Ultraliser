/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#include "BitArray.h"
#include "BitArray.hh"
#include <utilities/Utilities.h>

namespace Ultraliser
{

BitArray::BitArray(const uint64_t numBits)
    : _numBits(numBits)
{
    // Allocate space for bit array
    _numBytes = BITS_TO_CHARS(numBits);
    _data = new uint8_t[_numBytes];

    // Set all bits to 0
    std::fill_n(_data, _numBytes, 0);
}

BitArray::BitArray(uint8_t *array, const uint64_t numBits)
    : _numBits(numBits)
{
    _numBytes = BITS_TO_CHARS(numBits);
    _data = new uint8_t[_numBytes];

    std::memcpy(_data, array, _numBytes);
}

void BitArray::print(std::ostream &outStream)
{
    uint64_t size = BITS_TO_CHARS(_numBits);

    outStream.width(2);
    outStream.fill('0');

    // First byte
    outStream << std::uppercase << std::hex << I2UI64(_data[0]);

    for (uint64_t i = 1; i < size; i++)
    {
        // Remaining bytes with a leading space
        outStream << " ";
        outStream.width(2);
        outStream.fill('0');
        outStream << (_data[i]);
    }

    outStream << std::dec;
}

void BitArray::writeBIN(const std::string &prefix)
{
    std::ofstream outStream;
    outStream.open(prefix.c_str());

    uint64_t size = BITS_TO_CHARS(_numBits);

    for (uint64_t i = 0; i < size; i++)
        outStream << (_data[i]);

    outStream.close();
}

void BitArray::setAll()
{
    uint64_t size = BITS_TO_CHARS(_numBits);

    // Set bits in all bytes to 1
    std::fill_n(_data, size, UCHAR_MAX);

    // Zero any spare bits so increment and decrement are consistent
    uint64_t bits = _numBits % CHAR_BIT;
    if (bits != 0)
    {
        uint8_t mask = I2UI8(UCHAR_MAX << (CHAR_BIT - bits));
        _data[BIT_CHAR(_numBits - 1)] = mask;
    }
}

void BitArray::clearAll()
{
    uint64_t size = BITS_TO_CHARS(_numBits);

    // Set bits in all bytes to 0
    std::fill_n(_data, size, 0);
}

void BitArray::setBit(const uint64_t bit)
{
    if (_numBits <= bit)
        return; // Bit out of range

    _data[BIT_CHAR(bit)] |= BIT_IN_CHAR(bit);
}

void BitArray::clearBit(const uint64_t bit)
{
    // Bit out of range
    if (_numBits <= bit)
        return;

    // Create a mask to zero out desired bit
    uint8_t mask = I2UI8(BIT_IN_CHAR(bit));
    mask = ~mask;

    _data[BIT_CHAR(bit)] &= mask;
}

BitArrayIndex BitArray::operator()(const uint64_t bit)
{
    BitArrayIndex result(this, bit);
    return result;
}

bool BitArray::operator[](const uint64_t bit) const
{
    return((_data[BIT_CHAR(bit)] & BIT_IN_CHAR(bit)) != 0);
}

bool BitArray::bit(const uint64_t index) const
{
    return((_data[BIT_CHAR(index)] & BIT_IN_CHAR(index)) != 0);
}

uint8_t BitArray::getByte(const uint64_t index) const
{
    return _data[index];
}

void BitArray::addByte(const uint64_t index, const uint8_t byte)
{
    _data[index] |= byte;
}

bool BitArray::operator==(const BitArray &other) const
{
    if (_numBits != other._numBits) return false;

    return (this->_data == other._data);
}

bool BitArray::operator!=(const BitArray &other) const
{
    if (_numBits != other._numBits) return true;

    return (this->_data != other._data);
}

bool BitArray::operator<(const BitArray &other) const
{
    if (_numBits != other._numBits) return false;

    return (this->_data < other._data);
}

bool BitArray::operator<=(const BitArray &other) const
{
    if (_numBits != other._numBits) return false;

    return (this->_data <= other._data);
}

bool BitArray::operator>(const BitArray &other) const
{
    if (_numBits != other._numBits) return false;

    return (this->_data > other._data);
}

bool BitArray::operator>=(const BitArray &other) const
{
    if (_numBits != other._numBits) return false;

    return (this->_data >= other._data);
}

BitArray BitArray::operator~() const
{
    BitArray result(this->_numBits);
    result = *this;
    result.Not();

    return result;
}

BitArray BitArray::operator&(const BitArray &other) const
{
    BitArray result(this->_numBits);
    result = *this;
    result &= other;

    return result;
}

BitArray BitArray::operator^(const BitArray &other) const
{
    BitArray result(this->_numBits);
    result = *this;
    result ^= other;

    return result;
}

BitArray BitArray::operator|(const BitArray &other) const
{
    BitArray result(this->_numBits);
    result = *this;
    result |= other;

    return result;
}

BitArray BitArray::operator<<(const uint64_t count) const
{
    BitArray result(this->_numBits);
    result = *this;
    result <<= count;

    return result;
}

BitArray BitArray::operator>>(const uint64_t count) const
{
    BitArray result(this->_numBits);
    result = *this;
    result >>= count;

    return result;
}

BitArray& BitArray::operator++()
{
    // Nothing to increment
    if (_numBits == 0)
        return *this;

    // Handle arrays that don't use every bit in the last character
    int64_t i = (_numBits % CHAR_BIT);    

    // Maximum value for current char
    uint8_t maxValue;

    // Least significant bit in current char
    uint8_t one;
    if (i != 0)
    {
        maxValue = I2UI8(UCHAR_MAX << (CHAR_BIT - i));
        one = I2UI8(1 << (CHAR_BIT - i));
    }
    else
    {
        maxValue = UCHAR_MAX;
        one = 1;
    }

    for (i = BIT_CHAR(_numBits - 1); i >= 0; i--)
    {
        if (_data[i] != maxValue)
        {
            _data[i] = _data[i] + one;
            return *this;
        }
        else
        {
            // Need to carry to next byte
            _data[i] = 0;

            // Remaining characters must use all bits
            maxValue = UCHAR_MAX;
            one = 1;
        }
    }

    return *this;
}

BitArray& BitArray::operator++(int)
{
    ++(*this);
    return *this;
}

BitArray& BitArray::operator--()
{
    // Nothing to decrement
    if (_numBits == 0)
        return *this;

    // Handle arrays that don't use every bit in the last character
    int64_t i = (_numBits % CHAR_BIT);
    uint8_t maxValue;     // Maximum value for current char
    uint8_t one;          // Least significant bit in current char
    if (i != 0)
    {
        maxValue = I2UI8(UCHAR_MAX << (CHAR_BIT - i));
        one = I2UI8(1 << (CHAR_BIT - i));
    }
    else
    {
        maxValue = UCHAR_MAX;
        one = 1;
    }

    for (i = BIT_CHAR(_numBits - 1); i >= 0; i--)
    {
        if (_data[i] >= one)
        {
            _data[i] = _data[i] - one;
            return *this;
        }
        else
        {
            // Need to borrow from the next byte
            _data[i] = maxValue;

            // Remaining characters must use all bits
            maxValue = UCHAR_MAX;
            one = 1;
        }
    }

    return *this;
}

BitArray& BitArray::operator--(int)
{
    --(*this);
    return *this;
}

BitArray& BitArray::operator=(const BitArray &src)
{
    // Don't do anything for a self assignment
    if (*this == src)
        return *this;

    // Don't do assignment with different array sizes
    if (_numBits != src._numBits)
        return *this;

    // Don't do assignment with unallocated array
    if ((_numBits == 0) || (src._numBits == 0))
        return *this;

    // Copy bits from source
    std::copy(src._data, &src._data[BITS_TO_CHARS(_numBits)], this->_data);
    return *this;
}

BitArray& BitArray::operator&=(const BitArray &src)
{
    // Don't do assignment with different array sizes
    if (_numBits != src._numBits)
        return *this;

    // AND array one uint8_t at a time
    for (uint64_t i = 0; i < BITS_TO_CHARS(_numBits); i++)
        _data[i] = _data[i] & src._data[i];

    return *this;
}

BitArray& BitArray::operator^=(const BitArray &src)
{
    // Don't do assignment with different array sizes
    if (_numBits != src._numBits)
        return *this;

    // XOR array one uint8_t at a time
    for (uint64_t i = 0; i < BITS_TO_CHARS(_numBits); i++)
        _data[i] = _data[i] ^ src._data[i];

    return *this;
}

BitArray& BitArray::operator|=(const BitArray &src)
{
    uint64_t size = BITS_TO_CHARS(_numBits);

    // Don't do assignment with different array sizes
    if (_numBits != src._numBits)
        return *this;

    // OR array one uint8_t at a time
    for (uint64_t i = 0; i < size; i++)
        _data[i] = _data[i] | src._data[i];

    return *this;
}

BitArray& BitArray::Not()
{
    // Don't do not with unallocated array
    if (_numBits == 0)
        return *this;

    // NOT array one uint8_t at a time
    for (uint64_t i = 0; i < BITS_TO_CHARS(_numBits); i++)
        _data[i] = ~_data[i];

    // Zero any spare bits so increment and decrement are consistent
    uint64_t bits = _numBits % CHAR_BIT;
    if (bits != 0)
    {
        uint8_t mask = I2UI8(UCHAR_MAX << (CHAR_BIT - bits));
        _data[BIT_CHAR(_numBits - 1)] &= mask;
    }

    return *this;
}

BitArray& BitArray::operator<<=(const uint64_t shifts)
{
    int chars = I2I32(shifts / CHAR_BIT); // Number of whole byte shifts

    if (shifts >= _numBits)
    {
        // All bits have been shifted off
        this->clearAll();
        return *this;
    }

    // First handle big jumps of bytes
    if (chars > 0)
    {
        int size = I2I32(BITS_TO_CHARS(_numBits));

        for (int i = 0; (i + chars) < size; i++)
            _data[i] = _data[i + chars];

        // Now zero out new bytes on the right.
        for (int i = size; chars > 0; chars--)
            _data[i - chars] = 0;
    }

    // Now we have at most CHAR_BIT - 1 bit shifts across the whole array
    for (int i = 0; i < I2I32(shifts % CHAR_BIT); i++)
    {
        for (uint64_t j = 0; j < BIT_CHAR(_numBits - 1); j++)
        {
            _data[j] <<= 1;

            // Handle shifts across byte bounds.
            if (_data[j + 1] & MS_BIT)
                _data[j] |= 0x01;
        }

        _data[BIT_CHAR(_numBits - 1)] <<= 1;
    }

    return *this;
}

BitArray& BitArray::operator>>=(const uint64_t shifts)
{
    // Number of whole byte shifts
    int chars = I2I32(shifts / CHAR_BIT);

    if (shifts >= _numBits)
    {
        // All bits have been shifted off
        this->clearAll();
        return *this;
    }

    // First handle big jumps of bytes.
    if (chars > 0)
    {
        for (int i = I2I32(BIT_CHAR(_numBits - 1)); (i - chars) >= 0; i--)
            _data[i] = _data[i - chars];

        // Now zero out new bytes on the right
        for (; chars > 0; chars--)
            _data[chars - 1] = 0;
    }

    // Now we have at most CHAR_BIT - 1 bit shifts across the whole array
    for (int i = 0; i < I2I32(shifts % CHAR_BIT); i++)
    {
        for (uint64_t j = BIT_CHAR(_numBits - 1); j > 0; j--)
        {
            _data[j] >>= 1;

            // Handle shifts across byte bounds
            if (_data[j - 1] & 0x01)
                _data[j] |= MS_BIT;
        }

        _data[0] >>= 1;
    }

    /// Zero any spare bits that are shifted beyond the end of the bit array
    /// so that increment and decrement are consistent.
    int j = I2I32(_numBits % CHAR_BIT);
    if (j != 0)
    {
        uint8_t mask = I2UI8(UCHAR_MAX << (CHAR_BIT - j));
        _data[BIT_CHAR(_numBits - 1)] &= mask;
    }

    return *this;
}

BitArrayIndex::BitArrayIndex(BitArray *array, const uint64_t index)
{
    _bitArray = array;
    _index = index;
}

void BitArrayIndex::operator=(const bool src)
{
    // No array
    if (_bitArray == nullptr)
        return;

    // Index is out of bounds
    if (_bitArray->size() <= _index)
        return;

    if (src)
        _bitArray->setBit(_index);
    else
        _bitArray->clearBit(_index);
}

}
