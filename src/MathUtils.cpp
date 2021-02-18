/*
 * MathUtils.cpp
 *
 * Created on: Apr 12, 2010
 * Copyright: 2010, Bryan Holland
 *
 * Copyright notice:
 *
 * This file is part of 'OFR Analysis Tool' (OFR-AT).
 *
 * OFR-AT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OFR-AT is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OFR-AT.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ----------------------------------------------------------------------
 *
 * Description: see MathUtils.h
 *
 */

//project includes
#include "MathUtils.h"

//C++ includes
#include <cmath>

namespace math {

	//rounds off a double to the nearest integer
	int round(double r) {
		return (int) (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
	}

	//gives pi to 17 decimal places
	double pi() {
		return 3.14159265358979323;
	}
}
