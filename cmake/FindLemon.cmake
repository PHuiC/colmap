# Copyright (c), ETH Zurich and UNC Chapel Hill.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# Find package module for Lemon library.
#
# The following variables are set by this module:
#
#   LEMON_FOUND: TRUE if Lemon is found.
#   lemon: Imported target to link against.
#
# The following variables control the behavior of this module:
#
# LEMON_INCLUDE_DIR_HINTS: List of additional directories in which to
#                              search for Lemon includes.
# LEMON_LIBRARY_DIR_HINTS: List of additional directories in which to
#                              search for Lemon libraries.

set(LEMON_INCLUDE_DIR_HINTS "" CACHE PATH "Lemon include directory")
set(LEMON_LIBRARY_DIR_HINTS "" CACHE PATH "Lemon library directory")

unset(LEMON_FOUND)

find_package(lemon CONFIG QUIET)
if(TARGET lemon)
    set(LEMON_FOUND TRUE)
    message(STATUS "Found Lemon")
    message(STATUS "  Target : lemon")
else()
	list(APPEND LEMON_CHECK_INCLUDE_DIRS
        ${LEMON_INCLUDE_DIR_HINTS}
        /usr/include
        /usr/local/include
        /opt/include
        /opt/local/include
    )
	
	list(APPEND LEMON_CHECK_LIBRARY_DIRS
        ${LEMON_LIBRARY_DIR_HINTS}
        /usr/lib
        /usr/local/lib
        /opt/lib
        /opt/local/lib
    )
	
	find_path(LEMON_INCLUDE_DIRS
        NAMES
        lemon/list_graph.h
        PATHS
        ${LEMON_CHECK_INCLUDE_DIRS})
		
	find_library(LEMON_LIBRARIES
        NAMES
        lemon
        PATHS
        ${LEMON_CHECK_LIBRARY_DIRS})
		
	if(LEMON_INCLUDE_DIRS AND LEMON_LIBRARIES)
        set(LEMON_FOUND TRUE)
        message(STATUS "Found Lemon")
        message(STATUS "  Includes : ${LEMON_INCLUDE_DIRS}")
        message(STATUS "  Libraries : ${LEMON_LIBRARIES}")
	endif()
	
	add_library(lemon INTERFACE IMPORTED)
    target_include_directories(
		lemon INTERFACE ${LEMON_INCLUDE_DIRS})
    target_link_libraries(
		lemon INTERFACE ${LEMON_LIBRARIES})
endif()

if(NOT LEMON_FOUND AND LEMON_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find Lemon")
endif()