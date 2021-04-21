/// Example plugin for CRPropa.
///
/// Please consider sharing the awesome plugin with you fellow researchers by
/// creating a eperate repository for your project. We maintain a list of
/// plugins to CRPropa on our webpage and are happy to add a link to your
/// project, just send us: (name of the plugin, short description, url)

#include <radiopropa/Module.h>
#include <radiopropa/Source.h>


/// A custom C++ module
class MyModule : public radiopropa::Module
{
  public:
    /// The parent's constructor need to be called on initialization!
  	MyModule();
  	void process(radiopropa::Candidate *candidate) const;
};


/// A custom source feature
class AddMyProperty: public radiopropa::SourceFeature
{
  public:
    /// The parent's constructor need to be called on initialization!
  	AddMyProperty();
  	void prepareCandidate(radiopropa::Candidate &candidate) const;
};
