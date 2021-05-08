//////////////////////////
// parser.hpp
//////////////////////////
//
// Parser for settings file
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen
// Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef PARSER_HEADER
#define PARSER_HEADER

#include "metadata.hpp"

namespace gevolution
{

struct parameter
{
    char name[PARAM_MAX_LENGTH];
    char value[PARAM_MAX_LENGTH];
    bool used;
};

//////////////////////////
// readline
//////////////////////////
// Description:
//   reads a line of characters and checks if it declares a parameter; if yes,
//   i.e. the line has the format "<parameter name> = <parameter value>" (with
//   an optional comment added, preceded by a hash-symbol, '#'), the parameter
//   name and value are copied to the corresponding arrays, and 'true' is
//   returned. If the format is not recognized (or the line is commented using
//   the hash-symbol) 'false' is returned instead.
//
// Arguments:
//   line       string containing the line to be read
//   pname      will contain the name of the declared parameter (if found)
//   pvalue     will contain the value of the declared parameter (if found)
//
// Returns:
//   'true' if a parameter is declared in the line, 'false' otherwise.
//
//////////////////////////

bool readline (char *line, char *pname, char *pvalue);

//////////////////////////
// loadParameterFile
//////////////////////////
// Description:
//   loads a parameter file and creates an array of parameters declared therein
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     will contain the array of parameters (memory will be allocated)
//
// Returns:
//   number of parameters defined in the parameter file (= length of parameter
//   array)
//
//////////////////////////

int loadParameterFile (const char *filename, parameter *&params);

//////////////////////////
// saveParameterFile
//////////////////////////
// Description:
//   saves a parameter file
//
// Arguments:
//   filename   string containing the path to the parameter file
//   params     array of parameters
//   numparam   length of parameter array
//   used_only  if 'true', only the used parameters will be written (default)
//
// Returns:
//
//////////////////////////

void saveParameterFile (const char *filename, parameter *params,
                        const int numparam, bool used_only = true);

//////////////////////////
// parseParameter (int)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value
//   as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter
//   value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     int &pvalue);

//////////////////////////
// parseParameter (long)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value
//   as integer
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to integer which will contain the parsed parameter
//   value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     long &pvalue);

//////////////////////////
// parseParameter (double)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses its value
//   as double
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     reference to double which will contain the parsed parameter
//   value (if found)
//
// Returns:
//   'true' if parameter is found and parsed successfully, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     double &pvalue);

//////////////////////////
// parseParameter (char *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and retrieves its
//   value as string
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     character string which will contain a copy of the parameter
//   value (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     char *pvalue);

//////////////////////////
// parseParameter (double *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a
//   list of comma-separated double values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of double values which will contain the list of parsed
//   parameters (if found) nmax       maximum size of array; will be set to the
//   actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     double *pvalue, int &nmax);

//////////////////////////
// parseParameter (int *)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a
//   list of comma-separated integer values
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of integer values which will contain the list of parsed
//   parameters (if found) nmax       maximum size of array; will be set to the
//   actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     int *pvalue, int &nmax);

//////////////////////////
// parseParameter (char **)
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a
//   list of comma-separated strings
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     array of character strings which will contain the list of
//   parsed parameters (if found) nmax       maximum size of array; will be set
//   to the actual size at return
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseParameter (parameter *&params, const int numparam, const char *pname,
                     char **pvalue, int &nmax);

//////////////////////////
// parseFieldSpecifiers
//////////////////////////
// Description:
//   searches parameter array for specified parameter name and parses it as a
//   list of comma-separated field specifiers
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   pname      name of parameter to search for
//   pvalue     integer which will contain the binary-encoded list of parsed
//   specifiers (if found)
//
// Returns:
//   'true' if parameter is found, 'false' otherwise
//
//////////////////////////

bool parseFieldSpecifiers (parameter *&params, const int numparam,
                           const char *pname, int &pvalue);

//////////////////////////
// parseMetadata
//////////////////////////
// Description:
//   parses all metadata from the parameter array
//
// Arguments:
//   params     array of parameters
//   numparam   length of parameter array
//   sim        reference to metadata stucture (holds simulation parameters)
//   cosmo      reference to cosmology structure (holds cosmological
//   parameters) ic         reference to icsettings structure (holds settings
//   for IC generation)
//
// Returns:
//   number of parameters parsed
//
//////////////////////////

#ifndef LATFIELD2_HPP
#define COUT std::cout
#endif

int parseMetadata (parameter *&params, const int numparam, metadata &sim,
                   cosmology &cosmo, icsettings &ic);
}
#endif
