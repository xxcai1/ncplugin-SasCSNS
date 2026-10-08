#ifndef NCPlugin_Factory_hh
#define NCPlugin_Factory_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //Factory which implements logic of how the physics model provided by the
  //plugin should be combined with existing models in NCrystal:

  class PluginFactory final : public NC::FactImpl::ScatterFactory {
  public:
    const char * name() const noexcept override;
    NC::Priority query( const NC::FactImpl::ScatterRequest& ) const override;
    NC::ProcImpl::ProcPtr produce( const NC::FactImpl::ScatterRequest& ) const override;
    //The SANS model needs the combined (multi-phase) Info object in order to
    //compute the scattering length density of the involved phases:
    NC::FactImpl::MultiPhaseCapability multiPhaseCapability() const override
    { return NC::FactImpl::MultiPhaseCapability::Both; }
  };
}

#endif
