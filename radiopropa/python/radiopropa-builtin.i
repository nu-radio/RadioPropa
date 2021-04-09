/* CRPropa3 SWIG interface (for Python) */

/* Content:
 *
 * 1. SWIG settings and workarounds
 * 2. SWIG and CRPropa headers
 * 3. Pretty print for Python
 * 4. Magnetic Lens and Particle Maps Container
 *
 */


/* 1. SWIG settings and workarounds */
%include "1_swig.i"

/* 2: SWIG and CRPropa headers */
%include "2_headers.i"

/* 3. Pretty print for Python */

%define __REPR__( classname )
%feature("python:slot", "tp_str", functype="reprfunc") classname::repr();
%feature("python:slot", "tp_repr", functype="reprfunc") classname::repr();

%extend classname {
        const std::string repr() {
            return $self->getDescription();
        }
}
%enddef

%define VECTOR3__REPR__( classname )
%feature("python:slot", "tp_str", functype="reprfunc") classname::repr();
%feature("python:slot", "tp_repr", functype="reprfunc") classname::repr();

%exception classname::__getitem__ {
  try {
        $action
  }
  catch (RangeError) {
        SWIG_exception(SWIG_IndexError, "Index out of bounds");
        return NULL;
  }

}

/*
%extend classname {
        const std::string repr() {
            char buffer[1024];
            sprintf( buffer, "Vector(%.6G, %.6G, %.6G)", $self->x, $self->y, $self->z );
            return buffer;
        }
        double __getitem__(size_t i) {
                if(i == 0)
                        return $self->getX();
                if(i == 1)
                        return $self->getY();
                if(i == 2)
                        return $self->getZ();
                throw RangeError();
        }
}
*/


%template(Vector3d) radiopropa::Vector3<double>;
%template(Vector3f) radiopropa::Vector3<float>;

%enddef

/* Division of vector fix #34 */
%feature("python:slot", "nb_divide", functype="binaryfunc") *::operator/;

%include "3_repr.i"
