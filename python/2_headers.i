/* 2: CRPropa headers and Python extensions */

/* Python slots */
%feature("python:slot", "sq_length", functype="lenfunc") __len__;
%feature("python:slot", "mp_subscript", functype="binaryfunc") __getitem__;
%feature("python:slot", "tp_iter", functype="unaryfunc") __iter__;
#ifdef SWIG_PYTHON3
%feature("python:slot", "tp_iternext", functype="iternextfunc") __next__;
#else
%feature("python:slot", "tp_iternext", functype="iternextfunc") next;
#endif

/* Include headers */

#ifdef CRPROPA_HAVE_QUIMBY
%import (module="quimby") "quimby/Referenced.h"
%import (module="quimby") "quimby/Vector3.h"
//%import (module="quimby") quimby.i
#endif

#ifdef CRPROPA_HAVE_SAGA
%import (module="saga") saga.i
#endif

%{
#include "RadioPropa.h"
%}

%{
#include <iostream>
#include <iomanip>
%}

%ignore operator<<;
%ignore operator>>;
%ignore *::operator=;
%ignore operator radiopropa::Source*;
%ignore operator radiopropa::SourceList*;
%ignore operator radiopropa::SourceInterface*;
%ignore operator radiopropa::SourceFeature*;
%ignore operator radiopropa::Candidate*;
%ignore operator radiopropa::Module*;
%ignore operator radiopropa::ModuleList*;
%ignore operator radiopropa::ScalarField*;
%ignore operator radiopropa::Observer*;
%ignore operator radiopropa::ObserverFeature*;
%ignore operator radiopropa::ParticleCollector*;
%ignore radiopropa::TextOutput::load;

%feature("ref")   radiopropa::Referenced "$this->addReference();"
%feature("unref") radiopropa::Referenced "$this->removeReference();"


%include "radiopropa/Logging.h"
%include "radiopropa/Vector3.h"
%include "radiopropa/Referenced.h"
%include "radiopropa/Units.h"
%include "radiopropa/Common.h"
%include "radiopropa/Cosmology.h"
%include "radiopropa/Random.h"
%include "radiopropa/ParticleState.h"
%include "radiopropa/Geometry.h"
%include "radiopropa/Version.h"

%import "radiopropa/Variant.h"

/* override Candidate::getProperty() */
%ignore radiopropa::Candidate::getProperty(const std::string &) const;

%nothread; /* disable threading for extend*/
%extend radiopropa::Candidate {
    PyObject * getProperty(PyObject * name){

        std::string input;

        if (PyUnicode_Check(name)){
          #ifdef SWIG_PYTHON3
          // test on PY_MAJOR_VERSION >= 3 wont work with swig
              input = PyUnicode_AsUTF8(name);
          #else
              PyObject *s =  PyUnicode_AsUTF8String(name);
              input = PyString_AsString(s);
          #endif
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check(name)){
            input = PyString_AsString(name);
        }
        #endif
        else {
            std::cerr << "ERROR: The argument of getProperty() must be a string/unicode object!" << std::endl;
            return NULL;
        }

        radiopropa::Variant value = $self->getProperty(input);

        // implement this conversion here and not in the Variant as
        // __asPythonObject, as extensions cannot be called from extension.
        if (! value.isValid())
        {
          Py_INCREF(Py_None);
          return Py_None;
        }
        else if (value.getTypeInfo() == typeid(bool))
        {
         if(value.toBool())
         {
          Py_RETURN_TRUE;
         }
         else
         {
          Py_RETURN_FALSE;
         }
        }
        // convert all integer types to python long
        else if (value.getTypeInfo() == typeid(char))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(unsigned char))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int16_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint16_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int32_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint32_t))
        {
          return PyInt_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(int64_t))
        {
          return PyLong_FromLong(value.toInt64());
        }
        else if (value.getTypeInfo() == typeid(uint64_t))
        {
          return PyLong_FromUnsignedLong(value.toInt64());
        }
        // convert float and double to pyfloat which is double precision
        else if (value.getTypeInfo() == typeid(float))
        {
          return PyFloat_FromDouble(value.toDouble());
        }
        else if (value.getTypeInfo() == typeid(double))
        {
          return PyFloat_FromDouble(value.toDouble());
        }
        else if (value.getTypeInfo() == typeid(std::string))
        {
        #ifdef SWIG_PYTHON3
          return PyUnicode_FromString(value.toString().c_str());
        #else
          return PyString_FromString(value.toString().c_str());
        #endif
        }

        std::cerr << "ERROR: Unknown Type" << std::endl;
        return NULL;
    }


    PyObject * setProperty(PyObject * name, PyObject * value){

        std::string input;

        if (PyUnicode_Check(name)){
          #ifdef SWIG_PYTHON3
              input = PyUnicode_AsUTF8(name);
          #else
              input = PyUnicode_AS_DATA(name);
              PyObject *s =  PyUnicode_AsUTF8String(name);
              input = PyString_AsString(s);
          #endif
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( name )){
            input = PyString_AsString( name );
        }
        #endif
        else {
            std::cerr << "ERROR: The argument of setProperty() must be a string/unicode object!" << std::endl;
            return NULL;
        }


        if (value == Py_None)
        {
          $self->setProperty(input, radiopropa::Variant());
        Py_RETURN_TRUE;
        }
        else if (PyBool_Check(value))
        {
         if(value == Py_True)
         {
          $self->setProperty(input, true);
         }
         else
         {
          $self->setProperty(input, false);
         }
          Py_RETURN_TRUE;
        }
        else if (PyInt_Check(value))
        {
          $self->setProperty(input, radiopropa::Variant::fromInt32(PyInt_AsLong(value)));
          Py_RETURN_TRUE;
        }
        else if (PyLong_Check(value))
        {
          $self->setProperty(input, radiopropa::Variant::fromUInt64(PyLong_AsLong(value)));
          Py_RETURN_TRUE;
        }
        else if (PyFloat_Check(value))
        {
          $self->setProperty(input, radiopropa::Variant::fromDouble(PyFloat_AsDouble(value)));
          Py_RETURN_TRUE;
        }
        else if (PyUnicode_Check(value)){
        #ifdef SWIG_PYTHON3
          $self->setProperty(input, PyUnicode_AsUTF8(value));
        #else
          PyObject *s =  PyUnicode_AsUTF8String(value);
          $self->setProperty(input, PyString_AsString(s));
        #endif
          Py_RETURN_TRUE;
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( value))
        {
          $self->setProperty(input, PyString_AsString(value));
          Py_RETURN_TRUE;
        }
        #endif
        else
        {
          PyObject *t = PyObject_Str(PyObject_Type(value));
          std::string ot;

          #ifdef SWIG_PYTHON3
            ot = PyUnicode_AsUTF8(t);
          #else
            ot = PyString_AsString(t);
          #endif
          std::cerr << "ERROR: Unknown Type: " << ot << std::endl;
          return NULL;
        }
    }
};
%thread; /* reenable threading */


%template(CandidateVector) std::vector< radiopropa::ref_ptr<radiopropa::Candidate> >;
%template(CandidateRefPtr) radiopropa::ref_ptr<radiopropa::Candidate>;
%include "radiopropa/Candidate.h"


%template(ModuleRefPtr) radiopropa::ref_ptr<radiopropa::Module>;
%template(stdModuleList) std::list< radiopropa::ref_ptr<radiopropa::Module> >;

%feature("director") radiopropa::Module;
%feature("director") radiopropa::AbstractCondition;
%include "radiopropa/Module.h"

%implicitconv radiopropa::ref_ptr<radiopropa::ScalarField>;
%template(ScalarFieldRefPtr) radiopropa::ref_ptr<radiopropa::ScalarField>;
%include "radiopropa/ScalarField.h"

%include "radiopropa/Grid.h"
%include "radiopropa/GridTools.h"

%implicitconv radiopropa::ref_ptr<radiopropa::Grid<radiopropa::Vector3<float> > >;
%template(VectorGridRefPtr) radiopropa::ref_ptr<radiopropa::Grid<radiopropa::Vector3<float> > >;
%template(VectorGrid) radiopropa::Grid<radiopropa::Vector3<float> >;

%implicitconv radiopropa::ref_ptr<radiopropa::Grid<float> >;
%template(ScalarGridRefPtr) radiopropa::ref_ptr<radiopropa::Grid<float> >;
%template(ScalarGrid) radiopropa::Grid<float>;

%include "radiopropa/EmissionMap.h"
%implicitconv radiopropa::ref_ptr<radiopropa::EmissionMap>;
%template(EmissionMapRefPtr) radiopropa::ref_ptr<radiopropa::EmissionMap>;
%implicitconv radiopropa::ref_ptr<radiopropa::CylindricalProjectionMap>;
%template(CylindricalProjectionMapRefPtr) radiopropa::ref_ptr<radiopropa::CylindricalProjectionMap>;

%include "radiopropa/module/BreakCondition.h"
%include "radiopropa/module/Boundary.h"

%feature("director") radiopropa::Observer;
%feature("director") radiopropa::ObserverFeature;
%include "radiopropa/module/Observer.h"
%include "radiopropa/module/Discontinuity.h"
%include "radiopropa/module/SimplePropagation.h"
%include "radiopropa/module/PropagationCK.h"

%ignore radiopropa::Output::enableProperty(const std::string &property, const Variant& defaultValue, const std::string &comment = "");
%extend radiopropa::Output{
  PyObject * enableProperty(const std::string &name, PyObject* defaultValue, const std::string &comment="")
  {

       if (defaultValue == Py_None)
        {
          Py_RETURN_TRUE;
        }
        else if (PyBool_Check(defaultValue))
        {
         if(defaultValue == Py_True)
         {
          $self->enableProperty(name, true, comment);
         }
         else
         {
          $self->enableProperty(name, false, comment);
         }
          Py_RETURN_TRUE;
        }
        else if (PyInt_Check(defaultValue))
        {
          $self->enableProperty(name, radiopropa::Variant::fromInt32(PyInt_AsLong(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyLong_Check(defaultValue))
        {
          $self->enableProperty(name, radiopropa::Variant::fromInt64(PyLong_AsLong(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyFloat_Check(defaultValue))
        {
          $self->enableProperty(name, radiopropa::Variant::fromDouble(PyFloat_AsDouble(defaultValue)), comment);
          Py_RETURN_TRUE;
        }
        else if (PyUnicode_Check(defaultValue)){
        #ifdef SWIG_PYTHON3
          std::string ss = PyUnicode_AsUTF8(defaultValue);
        #else
          PyObject *s =  PyUnicode_AsUTF8String(defaultValue);
          std::string ss = PyString_AsString(s);
        #endif
          $self->enableProperty(name, ss, comment);
          Py_RETURN_TRUE;
        }
        #ifndef SWIG_PYTHON3
        else if (PyString_Check( defaultValue))
        {
          std::string ss = PyString_AsString(defaultValue);
          $self->enableProperty(name, ss, comment);
          Py_RETURN_TRUE;
        }
        #endif
        else
        {
          PyObject *t = PyObject_Str(PyObject_Type(defaultValue));
          std::string ot;

          #ifdef SWIG_PYTHON3
            ot = PyUnicode_AsUTF8(t);
          #else
            ot = PyString_AsString(t);
          #endif
          std::cerr << "ERROR: Unknown Type: " << ot << std::endl;
          return NULL;
        }

  }
}


%include "radiopropa/module/Output.h"
%include "radiopropa/module/TextOutput.h"

%include "radiopropa/module/HDF5Output.h"
%include "radiopropa/module/OutputShell.h"
%include "radiopropa/module/OutputROOT.h"

%template(IntSet) std::set<int>;
%include "radiopropa/module/Tools.h"

%template(SourceInterfaceRefPtr) radiopropa::ref_ptr<radiopropa::SourceInterface>;
%feature("director") radiopropa::SourceInterface;
%template(SourceFeatureRefPtr) radiopropa::ref_ptr<radiopropa::SourceFeature>;
%feature("director") radiopropa::SourceFeature;
%include "radiopropa/Source.h"

%inline %{
class ModuleListIterator {
  public:
        ModuleListIterator(
                radiopropa::ModuleList::iterator _cur,
                radiopropa::ModuleList::iterator _end) : 
                        cur(_cur), end(_end) {}
        ModuleListIterator* __iter__() { return this; }
        radiopropa::ModuleList::iterator cur;
        radiopropa::ModuleList::iterator end;
  };
%}

%extend ModuleListIterator {
#ifdef SWIG_PYTHON3
  radiopropa::ref_ptr<radiopropa::Module>& __next__() {
#else
  radiopropa::ref_ptr<radiopropa::Module>& next() {
#endif
    if ($self->cur != $self->end) {
        return *$self->cur++;
    }
    throw StopIterator();
  }
}

%extend radiopropa::ModuleList {
  ModuleListIterator __iter__() {
        return ModuleListIterator($self->begin(), $self->end()); 
  }
  radiopropa::ref_ptr<radiopropa::Module> __getitem__(size_t i) {
        if (i >= $self->size()) {
                throw RangeError();
        }
        return (*($self))[i];
  }
  size_t __len__() {
        return $self->size();
  }
};

%template(ModuleListRefPtr) radiopropa::ref_ptr<radiopropa::ModuleList>;
%include "radiopropa/ModuleList.h"

%template(ParticleCollectorRefPtr) radiopropa::ref_ptr<radiopropa::ParticleCollector>;

%inline %{
class ParticleCollectorIterator {
  public:
        ParticleCollectorIterator(
                radiopropa::ParticleCollector::iterator _cur,
                radiopropa::ParticleCollector::iterator _end) : 
                        cur(_cur), end(_end) {}
        ParticleCollectorIterator* __iter__() { return this; }
        radiopropa::ParticleCollector::iterator cur;
        radiopropa::ParticleCollector::iterator end;
  };
%}

%extend ParticleCollectorIterator {
#ifdef SWIG_PYTHON3
  radiopropa::ref_ptr<radiopropa::Candidate>& __next__() {
#else
  radiopropa::ref_ptr<radiopropa::Candidate>& next() {
#endif
    if ($self->cur != $self->end) {
        return *$self->cur++;
    }
    throw StopIterator();
  }
}

%extend radiopropa::ParticleCollector {
  ParticleCollectorIterator __iter__() {
        return ParticleCollectorIterator($self->begin(), $self->end()); 
  }
  radiopropa::ref_ptr<radiopropa::Candidate> __getitem__(size_t i) {
        if (i >= $self->size()) {
                throw RangeError();
        }
        return (*($self))[i];
  }
  std::vector< radiopropa::ref_ptr<radiopropa::Candidate> > __getitem__(PyObject *param) {
        std::vector< radiopropa::ref_ptr<radiopropa::Candidate> > result;

        if (PySlice_Check(param)) {
                Py_ssize_t len = 0, start = 0, stop = 0, step = 0, slicelength = 0, i = 0;
                len = $self->size();

                #ifdef SWIG_PYTHON3
                    PySlice_GetIndicesEx(param, len, &start, &stop, &step, &slicelength);
                #else
                    PySlice_GetIndicesEx((PySliceObject*)param, len, &start, &stop, &step, &slicelength);
                #endif
                
                for(radiopropa::ParticleCollector::iterator itr = $self->begin(); itr != $self->end(); ++itr){
                        if( i >= start && i < stop){
                                result.push_back(itr->get());
                        }
                        ++i;
                }
                return result;
        } else {
                throw RangeError();
        }        
  }
  size_t __len__() {
        return $self->size();
  }
};

%include "radiopropa/module/ParticleCollector.h"
