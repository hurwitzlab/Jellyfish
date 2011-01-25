/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JELLYFISH_FLOATS_HPP__
#define __JELLYFISH_FLOATS_HPP__

#include <stdint.h>

namespace jellyfish {
  class Float {
  public:
    typedef uint32_t bits_t;

  private:
    union float_int {
      float      fv;
      bits_t     iv;
      float_int(float v) : fv(v) {}
      float_int(bits_t v) : iv(v) {}
    };
    float_int v;

  public:
    Float() : v(0.0f) {}
    Float(float _v) : v(_v) {}
    Float(bits_t _v) : v(_v) {}

    static const Float zero;
    static const Float one;

    const Float operator+(const Float y) const {
      return Float(v.fv + y.v.fv);
    }

    bits_t bits() const { return v.iv; };
    float to_float() const { return v.fv; };
  };
}

#endif