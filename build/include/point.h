/**
* @file   point.h
* @author   T.Sato
* @date   2015.04.28
* @version   1.0
*/

#ifndef __POINT__
#define __POINT__

/**
* @struct point_t
* @brief this structure has propaties of a point
*/
template<typename data_t>
struct point_t
{

	size_t index;			//!< index of data
	data_t* addressOfpoint;	//!< head address of a data
	double distance;		//!< ditance from query

	/**
	* @brief default constructor
	*/
	point_t()
	{}

	/**
	* @brief constructor
	*/
	point_t(
		const size_t& index,	//!< [in] the index of point 
		data_t* addressOfpoint, //!< [in] the address of point 
		const double& distance	//!< [in] the distance from query
		)
		: index(index)
		, addressOfpoint(addressOfpoint)
		, distance(distance)
	{}

	/**
	* @brief constructor
	*/
	point_t(
		const size_t& index,	//!< [in] the index of point 
		const double& distance	//!< [in] the distance from query
		)
		: index(index)
		, distance(distance)
	{}

	/**
	* @brief constructor
	*/
	point_t(
		const size_t& index,	//!< [in] the index of point 
		data_t* addressOfpoint	//!< [in] the address of point
		)
		: index(index)
		, addressOfpoint(addressOfpoint)
	{}

	/**
	* @brief set member variables
	*/
	void setMemberVariable(
		const size_t& index,	//!< [in] the index of point
		data_t* addressOfpoint,	//!< [in] the address of point
		const double& distance	//!< [in] the distance from query
		)
	{
		this->index = index;
		this->addressOfpoint = addressOfpoint;
		this->distance = distance;
	}

	/**
	* @brief compare the distance
	* @return is my distance lessor than e's ?
	*/
	bool operator <(
		const point_t &e	//!< [in] compare object
		) const
	{
		return distance < e.distance;
	}

	/**
	* @brief compare the distance
	* @return is my distance equal to e's ?
	*/
	bool operator ==(
		const point_t &e	//!< [in] compare object
		) const
	{
		return index == e.index;
	}
};

#endif
