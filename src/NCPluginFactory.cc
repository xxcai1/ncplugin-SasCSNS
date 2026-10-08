
#include "NCPluginFactory.hh"
#include "NCDirectLoad2D.hh"
#include "NCSansIsotropic.hh"
#include "NCSansModelPicker.hh"

namespace NCPluginNamespace {

  class PluginScatter final : public NC::ProcImpl::ScatterIsotropicMat {
  public:

    //The factory wraps our custom PhysicsModel helper class in an NCrystal API
    //Scatter class.

    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "Model"; }
    PluginScatter( SansIsotropic && pm ) : m_pm(std::move(pm)) {}

    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin) const override
    {
      return NC::CrossSect{ m_pm.calcCrossSection(ekin.dbl()) };
    }

    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG& rng, NC::NeutronEnergy ekin ) const override
    {
      auto outcome = m_pm.sampleScatteringEvent( rng, ekin.dbl() );
      return { NC::NeutronEnergy{outcome.ekin_final}, NC::CosineScatAngle{outcome.mu} };
    }

  private:
    SansIsotropic m_pm;
  };

}

const char * NCP::PluginFactory::name() const noexcept
{
  //Factory name. Keep this standardised form please:
  return NCPLUGIN_NAME_CSTR "Factory";
}

NC::Priority NCP::PluginFactory::query( const NC::FactImpl::ScatterRequest& req ) const
{
  //Respect the "sans" parameter, which allows users to disable SANS models:
  if ( ! req.get_sans() )
    return NC::Priority::Unable;
  //The model below can only deal with single-phase Info objects (it reads its
  //@CUSTOM_SASCSNS section via Info::countCustomSections, which throws on
  //multi-phase Info objects). Decline multi-phase requests here so that such
  //materials (e.g. files with @OTHERPHASES, including those using the builtin
  //@CUSTOM_HARDSPHERESANS section) keep working, served by NCrystal's own
  //multi-phase capable SANS factories instead:
  if ( req.isMultiPhase() )
    return NC::Priority::Unable;
  if ( ! SansModelPicker::isApplicable(req.info()) )
    return NC::Priority::Unable;
  return NC::Priority{999};
}

NC::ProcImpl::ProcPtr NCP::PluginFactory::produce( const NC::FactImpl::ScatterRequest& req ) const
{
  //Tabulated anisotropic 2D SANS (DirectLoad2D): the material becomes
  //anisotropic and transport must ask for direction-dependent cross sections:
  if ( DirectLoad2D::isDirectLoad2D(req.info()) ) {
    return NC::makeSO<DirectLoad2DScatter>( DirectLoad2D::createFromInfo(req.info()) );
  }

  auto sc_ourmodel = NC::makeSO<PluginScatter>(SansModelPicker::createFromInfo(req.infoPtr()));
  return sc_ourmodel;

  //To add our SANS contribution on top of the standard NCrystal models, the
  //following could be used instead:
  // auto sc_std = globalCreateScatter( req );
  // return combineProcs( sc_std, sc_ourmodel );
}
