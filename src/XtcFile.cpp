/*
 * XtcFile.cpp
 *
 *  Created on: May 3, 2013
 *      Author: bholland
 */

#include "XtcFile.h"

#include <iostream>
#include <algorithm>
#include <climits>
#include <QtCore/QtDebug>

using namespace std;

//Need the magic number used by XTC files
const int XTC_MAGIC = 1995;
const int ATOMIC_DEGREES_OF_FREEDOM = 3;
const int BYTES_PER_FLOAT = 4;
const int MAX_SMALLER_VALUE = 0xffffff;
const int BITS_PER_QINT32 = 32;
const int BITS_PER_BYTE = 8;
const int LEAST_SIGBYTE_MASK = 0xff;
const uchar CHAR_SIGN_CONVERSION_BYTE_MASK = 0xff;
const int BINARY_ORDER_OF_MAGNITUDE = 2;
const double ANGSTROMS_PER_NM = 10;

const int MAXABS = INT_MAX - 2;

//Default constructor
XtcFile::XtcFile() {

	initializeMembers();
}

//Constructor that sets the QFile to the given file name, and opens a QDataStream
//using that QFile.
XtcFile::XtcFile(const QString& file_name, System* system) : QFile(file_name) {

	this->system = system;
	initializeMembers();
}

//getters and setters
uint XtcFile::numberOfAtoms() const {return num_atoms;}
float XtcFile::precision() const {return _precision;}
uint XtcFile::firstStep() const {return first_step;}
float XtcFile::startTime() const {return start_time;}
double XtcFile::timeStep() const {return time_step;}
uint XtcFile::currentStep() const {return cur_step;}
float XtcFile::currentTime() const {return cur_time;}

void XtcFile::setSystem(System* system) {

	this->system = system;
}

//Opens the file and the data stream using the current object. Don't want to do in
//constructor so a 'success' bool can be returned
bool XtcFile::openXtc() {

	if(exists()) {

		if(not open(QIODevice::ReadOnly))
			return false;

		stream.setDevice(this);
		stream.setFloatingPointPrecision(QDataStream::SinglePrecision); //xtc only reads in 4 byte floats
		return true;
	}

	return false;
}

//Reads in all of the header information from the XtcFile.  File must be open, or outputs
//an error message and does nothing.
bool XtcFile::readHeader() {

	//error checks
	if(not isOpen()) {

		errorStream("XTC file must be open to read header!");
		return false;

	} else if(atEnd()) {

		errorStream("Already at end of file, no header here!");
		return false;
	}

	//some dummy variables
	float dummy_float;
	int dummy_int;

	//make sure file is at beginning
	resetToStart();

	//check the "magic number", these things are stupid
	int magic;
	stream >> magic;

	if(magic != XTC_MAGIC) {

		errorStream("Wrong magic number, maybe not an XTC file or file corrupt?");
		return false;
	}

	//Read in some system info
	stream >> num_atoms;
	stream >> first_step;
	stream >> start_time;

	// skip the simulation box for the header;
	for(uint i = 0; i < 9; ++i)
		stream >> dummy_float;

	//Number of coordinate in the frame - should be the same as the number of atoms, redundant
	stream >> dummy_int;

	if(num_atoms <= 9){

		_precision = -1;					// no compression;
		coord_block_size = num_atoms * ATOMIC_DEGREES_OF_FREEDOM * BYTES_PER_FLOAT;
	}
	else {

		stream >> _precision;

		// skip the ranges and the small value (?)
		for(int i = 0; i < 7; ++i)
			stream >> dummy_int;

		stream >> coord_block_size;

		//adjust for padding to 4N XDR block size;
		coord_block_size = ((coord_block_size + 3) / 4) * BYTES_PER_FLOAT;
	}

	//reset to the beginning of the first frame
	resetToStart();
	return true;
}

//Reads the next frame from the current position of the file and sends the results out in a Qt signal.
bool XtcFile::readNextFrame() {

	//error checks - if error exists return empty XtcFrame
	if(not isOpen()) {

		errorStream("XTC file must be open to read frames!");
		return false;

	} else if(atEnd()) {

		errorStream("Reached EOF, frame returned is empty");
		return false;
	}

	//check the "magic number"
	int magic;
	stream >> magic;

	if(magic != XTC_MAGIC) {

		QString msg;
		msg.append("Wrong magic number at frame ");
		msg.append(QString::number(cur_step+1));
		msg.append(" , file might be corrupt!");
		errorStream(msg);
		return false;
	}

	int natoms;
	float v1, v2, v3;

	stream >> natoms;
	stream >> cur_step;
	stream >> cur_time;

	//calculate the time step if hasn't been done yet -> this should clearly be done on the second frame read
	if(not time_step_calculated and cur_time > start_time) {

		time_step = cur_time - start_time;
		time_step_calculated = true;
	}

	//space vectors - convert to angstroms, this is the standard unit of length in SIMtoEXP
	stream >> v1;
	stream >> v2;
	stream >> v3;
	Vector3D a(v1, v2, v3);
	a.scale(ANGSTROMS_PER_NM);

	stream >> v1;
	stream >> v2;
	stream >> v3;
	Vector3D b(v1, v2, v3);
	b.scale(ANGSTROMS_PER_NM);

	stream >> v1;
	stream >> v2;
	stream >> v3;
	Vector3D c(v1, v2, v3);
	c.scale(ANGSTROMS_PER_NM);

	//pass along the box info to the system
	system->updateMDBox(Matrix3D (a,b,c));

	//read all the coordinate for this frame and store them in the 'Frame' object
	uint num_read = readCompressedCoords();
	if(num_read == 0)
		return false;

	return true;
}

//private function that does the work of xdr3dfcoord from the Gromacs xdr package (?).  Reads
//in and uncompresses the coordinates. Returns the number of coordinates read
uint XtcFile::readCompressedCoords() {

	/*____________________________________________________________________________
	 |
	 | xdr3dfcoord - read or write compressed 3d coordinates to xdr file.
	 |
	 | this routine reads or writes (depending on how you opened the file with
	 | xdropen() ) a large number of 3d coordinates (stored in *fp).
	 | The number of coordinates triplets to write is given by *size. On
	 | read this number may be zero, in which case it reads as many as were written
	 | or it may specify the number if triplets to read (which should match the
	 | number written).
	 | Compression is achieved by first converting all floating numbers to integer
	 | using multiplication by *precision and rounding to the nearest integer.
	 | Then the minimum and maximum value are calculated to determine the range.
	 | The limited range of integers so found, is used to compress the coordinates.
	 | In addition the differences between succesive coordinates is calculated.
	 | If the difference happens to be 'small' then only the difference is saved,
	 | compressing the data even more. The notion of 'small' is changed dynamically
	 | and is enlarged or reduced whenever needed or possible.
	 | Extra compression is achieved in the case of GROMOS and coordinates of
	 | water molecules. GROMOS first writes out the Oxygen position, followed by
	 | the two hydrogens. In order to make the differences smaller (and thereby
	 | compression the data better) the order is changed into first one hydrogen
	 | then the oxygen, followed by the other hydrogen. This is rather special, but
	 | it shouldn't harm in the general case.
	 |
	 */

//	int xdr3dfcoord(FILE *fptr, valarray<float>& xx, size_t& size, float& precision)

	//get the number of coordinates and allocate the space in the vector to avoid numerous resizes
	int num_coords;
	stream >> num_coords;

	//deal with possible error
	if(num_coords < 0) {

		errorStream("Number of coordinates for frame is < 0, something is horribly wrong!");
		return 0;
	}

	//if there are < 10 particles, no compression
	float v1, v2, v3;
	if (num_coords <= 9) {

		//precision is meaningless without compression
		_precision = -1;

		//go through and add all the coordinates to the frame object
		for(int i = 0; i < num_coords; ++i) {

			stream >> v1;
			stream >> v2;
			stream >> v3;

			setCoords(v1, v2, v3);
		}

		return num_coords;
	}

	//otherwise, need to deal with compression
	stream >> _precision;

	int min_ix, min_iy, min_iz, max_ix, max_iy, max_iz;
	stream >> min_ix;
	stream >> min_iy;
	stream >> min_iz;

	stream >> max_ix;
	stream >> max_iy;
	stream >> max_iz;

	QVector<int> size_ints;
	size_ints << max_ix - min_ix + 1;
	size_ints << max_iy - min_iy + 1;
	size_ints << max_iz - min_iz + 1;

	//check if one of the sizes is too big to be multiplied
	quint32 bitsize_ix = 0, bitsize_iy = 0, bitsize_iz = 0;
	quint32 bitsize;

	if ((size_ints[0] bitor size_ints[1] bitor size_ints[2] ) > MAX_SMALLER_VALUE) {

		bitsize_ix = bitsizeQInt32(size_ints[0]);
		bitsize_iy = bitsizeQInt32(size_ints[1]);
		bitsize_iz = bitsizeQInt32(size_ints[2]);

		// Flag the use of large sizes
		bitsize = 0;
	}
	else {

		bitsize = bitsizeQInt32(size_ints);
	}

	int smallidx, smaller, smallnum, nbytes;

	//read the 'smallidx' value, not really sure what this is
	stream >> smallidx;

	smaller = magic_ints[max(FIRSTIDX, smallidx -1)] / 2;
	smallnum = magic_ints[smallidx] / 2;		// Factor for value/difference decision.

	QVector<int> sizesmall;
	sizesmall << magic_ints[smallidx] << magic_ints[smallidx] << magic_ints[smallidx];

	//Get the number of bytes to read in
	stream >> nbytes;

	// Read some bytes from the file - gets stored in member 'filed_data'
	readBytes(nbytes);

	// -----------------------------------------------------------------------------------------
	float inv_precision = 1.0 / _precision;

	int run = 0;						// ??? for water rearrangements?

	int cur_coord = 0; //counter
	//int fi = 0;	//float number counter - NOT USED BELOW

	//integer coordinates - TODO: vector note - originally 'size3' used to initialize, but
	//this doesn't seem to make sense as 'size3' was used to account for all 3 dimensions, whereas below
	//each dimension is initialized separately.  Use 'num_coords', but make sure to check!
	QVector<int> int_x(num_coords), int_y(num_coords), int_z(num_coords);

	int prev_int_x, prev_int_y, prev_int_z;  // previous integer coordinates;
	float float_x, float_y, float_z;

	quint32 count = 0, last_byte = 0, last_bit = 0;

	//Loop until all the vectors are read
	while (cur_coord < num_coords) {

		//use 'ic' to track 'integer array' position, starts at whatever the current coord is at
		int ic = cur_coord++;
		QVector<int> sizes;

		int ints_to_read = 3;

		// Decoding the binary stream we have just read as opaque data.
		//remember bitsize = 0 is for large values
		if (bitsize == 0) {

			// For larger values exceeding one integer space.
			int_x[ic] = receiveBits(bitsize_ix, count, last_byte, last_bit);
			int_y[ic] = receiveBits(bitsize_iy, count, last_byte, last_bit);
			int_z[ic] = receiveBits(bitsize_iz, count, last_byte, last_bit);
		}
		else {
			// For smaller values which fit into one integer value.
			receiveInts(ints_to_read, bitsize, size_ints, count, last_byte, last_bit);
			int_x[ic] = int_buffer[0];
			int_y[ic] = int_buffer[1];
			int_z[ic] = int_buffer[2];
		}

		int_x[ic] += min_ix;
		int_y[ic] += min_iy;
		int_z[ic] += min_iz;

		prev_int_x = int_x[ic];	// remember the 'previous' coords for the next iteration.
		prev_int_y = int_y[ic];
		prev_int_z = int_z[ic];

		quint32 flag_bitsize = 1;
		int flag = receiveBits(flag_bitsize, count, last_byte, last_bit);
		int is_smaller = 0;

		if(flag) {

			//if flag set, use next 5 bits to alter 'run' in some weird way.
			uint run_bitsize = 5;
			run = receiveBits(run_bitsize, count, last_byte, last_bit);
			is_smaller = run % 3;
			run -= is_smaller;
			is_smaller--;
		}

		if (run > 0) {

			ic++;

			for (int k = 0; k < run; k += 3) {

				receiveInts(ints_to_read, smallidx, sizesmall, count, last_byte, last_bit);
				int_x[ic] = int_buffer[0];
				int_y[ic] = int_buffer[1];
				int_z[ic] = int_buffer[2];

				cur_coord++;

				int_x[ic] += prev_int_x - smallnum;
				int_y[ic] += prev_int_y - smallnum;
				int_z[ic] += prev_int_z - smallnum;

				if (k == 0) {

					/* interchange first with second atom for better
					 * compression of water molecules
					 */
					swapValues(int_x[ic], prev_int_x);
					swapValues(int_y[ic], prev_int_y);
					swapValues(int_z[ic], prev_int_z);

					float_x = prev_int_x * inv_precision;
					float_y = prev_int_y * inv_precision;
					float_z = prev_int_z * inv_precision;

					setCoords(float_x, float_y, float_z);
				}
				else {

					prev_int_x = int_x[ic];
					prev_int_y = int_y[ic];
					prev_int_z = int_z[ic];
				}

			    float_x = int_x[ic] * inv_precision;
			    float_y = int_y[ic] * inv_precision;
			    float_z = int_z[ic] * inv_precision;

				setCoords(float_x, float_y, float_z);
			}
		}
		else {

		    float_x = int_x[ic] * inv_precision;
		    float_y = int_y[ic] * inv_precision;
		    float_z = int_z[ic] * inv_precision;

		    setCoords(float_x, float_y, float_z);
		}

		smallidx += is_smaller;

		if (is_smaller < 0) {

			smallnum = smaller;

			if (smallidx > FIRSTIDX) {

				smaller = magic_ints[smallidx - 1] /2;
			}
			else {

				smaller = 0;
			}
		}
		else {

			if (is_smaller > 0) {

				smaller = smallnum;
				smallnum = magic_ints[smallidx] / 2;
			}
		}

		sizesmall[0] = sizesmall[1] = sizesmall[2] = magic_ints[smallidx] ;

	} 							// end of the main while loop.
	// -----------------------------------------------------------------------------------------

	return num_coords;
}

//Moves the file cursor back to the beginning of the file
void XtcFile::resetToStart() {

	reset();
	cur_step = 0;
}

//Private function for determining the minimum number of bits required to store a qint32 type value
//Altered version of 'xdr' function from Gromacs 4.5
uint XtcFile::bitsizeQInt32(int value) const {

    int num = 1; //start with 31 zeros and 1
    int num_of_bits = 0;

    //count the bits needed
    while (value >= num and num_of_bits < BITS_PER_QINT32) {

    	++num_of_bits;
    	num <<= 1; //left-shift the single bit
    }

    return num_of_bits;
}

//Private function for determining minimum number of bits required to store an array of qint32 type values
//Altered version of 'xdr' function from Gromacs 4.5
uint XtcFile::bitsizeQInt32(const QVector<int>& sizes) const {

	const uint buffer = 32;
	QVector<int> bytes(buffer);
	quint32 num_of_bytes, num_of_bits, byte_cnt, tmp;

	//initialize number of bytes to 1 and bits to zero (so assume there will be at least one byte??)
	num_of_bytes = 1;
	bytes[0] = 1;
	num_of_bits = 0;

	int size = sizes.size();
	for(int i = 0; i < size; ++i) {

		tmp = 0;
		for(byte_cnt = 0; byte_cnt < num_of_bytes; ++byte_cnt) {

			tmp += bytes[byte_cnt] * sizes[i];
			bytes[byte_cnt] = tmp bitand LEAST_SIGBYTE_MASK;
			tmp >>= BITS_PER_BYTE;
		}

		while(tmp != 0) {

			bytes[byte_cnt++] = tmp bitand LEAST_SIGBYTE_MASK;
			tmp >>= BITS_PER_BYTE;
		}

		num_of_bytes = byte_cnt;
	}

	//remove the last byte counted and count the number of bits in it
	int ref_num = 1;
	--num_of_bytes;

	while (bytes[num_of_bytes] >= ref_num) {

		++num_of_bits;
		ref_num *= BINARY_ORDER_OF_MAGNITUDE;
	}

	return num_of_bits + num_of_bytes * BITS_PER_BYTE;
}

/*
 * decode number from buffer using specified number of bits
 * extract the number of bits from the vector and construct an integer from it. Return that value.
 *
 */

int XtcFile::receiveBits(quint32 num_of_bits, quint32& cnt, quint32& lastbyte, quint32& endbits) {

	//create a mask with all bits up to 'num_of_bits' set to 1
	int bit_mask = (1 << num_of_bits) - 1;
	int decoded_value = 0;

	//deal with all of the bytes first
	while (num_of_bits >= (quint32) BITS_PER_BYTE) {

		uchar cur_data = file_data[cnt++] bitand CHAR_SIGN_CONVERSION_BYTE_MASK;
		lastbyte = (lastbyte << BITS_PER_BYTE) bitor (cur_data); //increment the counter after the operation for next time
		num_of_bits -= BITS_PER_BYTE;
		decoded_value |=  (lastbyte >> endbits) << num_of_bits;
	}

	//now deal with the bits left over
	if (num_of_bits > 0) {

		//if the 'lastbits' don't cover the rest of bits needed, add a byte
		if (endbits < num_of_bits) {

			endbits += BITS_PER_BYTE;
//			qDebug() << "File data in if: " << (qint32) file_data[cnt];
			uchar cur_data = file_data[cnt++] bitand CHAR_SIGN_CONVERSION_BYTE_MASK;
			lastbyte = (lastbyte << BITS_PER_BYTE) bitor cur_data;
		}

		//now have enough, determine the new 'endbits' and the remaining number by shifting
		//'endbits' to the right and only keeping 'num_of_bits' through the mask
		endbits -= num_of_bits;
		decoded_value |= (lastbyte >> endbits) bitand bit_mask;
	}

	return (decoded_value bitand bit_mask);
}

/*
 * decode 'small' integers from the buffer array
 *
 * this routine is the inverse from sendints() and decodes the small integers
 * written to buffer by calculating the remainder and doing divisions with
 * the given sizes[]. You need to specify the total number of bits to be
 * used from the buffer in num_of_bits.
 *
 */

void XtcFile::receiveInts(quint32 num_of_ints, quint32 num_of_bits, QVector<int>& sizes,
		quint32& cnt, quint32& lastbyte, quint32& endbits) {

	//create an array for the 'bytes'
	int buffer_size = 32;
	QVector<int> bytes(buffer_size, 0);

	int num_of_bytes = 0;

	while(num_of_bits > (quint32) BITS_PER_BYTE) {

		bytes[num_of_bytes++] = receiveBits(BITS_PER_BYTE, cnt, lastbyte, endbits);
		num_of_bits -= BITS_PER_BYTE;
	}

	if (num_of_bits > 0) {

		bytes[num_of_bytes++] = receiveBits(num_of_bits, cnt, lastbyte, endbits);
	}

	for(quint32 i = num_of_ints-1; i > 0; --i) {

		int num = 0;

		for (int j = num_of_bytes-1; j >= 0; --j) {

			num = (num << BITS_PER_BYTE) bitor bytes[j];
			bytes[j] = num / sizes[i]; //integer division, so should truncate
			num -= bytes[j] * sizes[i];
		}
		int_buffer[i] = num;
	}

	int_buffer[0] = bytes[0] bitor (bytes[1] << BITS_PER_BYTE) bitor (bytes[2] << (2*BITS_PER_BYTE)) bitor (bytes[3] << (3*BITS_PER_BYTE));
}

QVector<char>& XtcFile::readBytes(int nbytes) {

	// It seems that XDR data can only be a multiple of 4 bytes,
	// so here we have to read some padding to 4n data length;
	int mod = (nbytes % 4);
	if(mod > 0)
		nbytes += (4 - mod);

	//need to preallocate since writing directly to the vector without it knowing about it (tricky!)
	//reserving is probably not good enough here though, as it will leave the size info out of sync
	if(file_data.size() < nbytes)
		file_data.resize(nbytes);

	stream.readRawData(&file_data[0], nbytes);

	return file_data;
}

//Private helper function for writing the given string to the standard err output
void XtcFile::errorStream(const QString& err) const {

	cerr << err.toStdString() << endl;
}

//Private function that add the magic integers to the container for reference later
void XtcFile::initializeMembers() {

	cur_step = 0;
	cur_time = 0;
	num_atoms = 0;
	_precision = 0;

	first_step = 0;
	start_time = 0;
	time_step = 0;
	time_step_calculated = false;

	coord_block_size = 0;

    magic_ints << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 <<
                  8 << 10 << 12 << 16 << 20 << 25 << 32 << 40 << 50 << 64 <<
                  80 << 101 << 128 << 161 << 203 << 256 << 322 << 406 << 512 << 645 <<
                  812 << 1024 << 1290 << 1625 << 2048 << 2580 << 3250 << 4096 << 5060 << 6501 <<
                  8192 << 10321 << 13003 << 16384 << 20642 << 26007 << 32768 << 41285 << 52015 << 65536 <<
                  82570 << 104031 << 131072 << 165140 << 208063 << 262144 << 330280 << 416127 << 524287 << 660561 <<
                  832255 << 1048576 << 1321122 << 1664510 << 2097152 << 2642245 << 3329021 << 4194304 << 5284491 << 6658042 <<
                  8388607 << 10568983 << 13316085 << 16777216;


    //find first index in magic number that is non-zero
    int i = 0;

    while(magic_ints[i] == 0 and i < magic_ints.size()) {
    	i++;
    }

    FIRSTIDX = i;
    LASTIDX = magic_ints.size()-1;

    int_buffer.resize(3);
}

//Private helper function that swaps two integer values
void XtcFile::swapValues(int& a, int& b) {

	int tmp = a;
	a = b;
	b = tmp;
}

//Private function that scales the coordinates to Angstroms and emits them to the ether
void XtcFile::setCoords(double x, double y, double z) {

	Vector3D coords(x, y, z);
	coords.scale(ANGSTROMS_PER_NM);
	system->updateNextCoord(coords);
}
