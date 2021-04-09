/* 3. Pretty print for Python */
/*  __repr__ << getDescription */

__REPR__( radiopropa::ParticleState );
__REPR__( radiopropa::Candidate );
__REPR__( radiopropa::Module );
__REPR__( radiopropa::ModuleList );
__REPR__( radiopropa::Source );
__REPR__( radiopropa::SourceList );
__REPR__( radiopropa::SourceFeature );
__REPR__( radiopropa::Observer );
__REPR__( radiopropa::ObserverFeature );

VECTOR3__REPR__( radiopropa::Vector3 );

%pythoncode %{
    DeflectionCK = PropagationCK  # legacy name
%}

%pythoncode %{
        __version__ = g_GIT_DESC 
%}

