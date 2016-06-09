/*
 * macros.hpp
 *
 *  Created on: Nov 13, 2014
 *      Author: rrsettgast
 */

#ifndef MACROS_HPP_
#define MACROS_HPP_

#define STRINGIZE_(x) #x
#define STRINGIZE(x) STRINGIZE_(x)

#define LOCATION __FILE__ ":" STRINGIZE(__LINE__)


#endif /* MACROS_HPP_ */
