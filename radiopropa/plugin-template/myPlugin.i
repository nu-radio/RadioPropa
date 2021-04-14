/* name of the plugin: myPlugin*/
%module(directors="1", threads="1", allprotected="1") myPlugin

/* Exceptions required */
%include "exception.i"

/*  define headers to include into the wrapper. These are the plugin headers
 *  and the RadioPropa headers.
 */
%{
#include "RadioPropa.h"
#include "myPlugin.h"
%}

/* import radiopropa in wrapper */
%import (module="radiopropa") "radiopropa.i"

/* include plugin parts to generate wrappers for */
%include "myPlugin.h"





