
#include "NCTestPlugin.hh"
#include "NCrystal/factories/NCFact.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

#include <cmath>
#include <string>

namespace NCPluginNamespace {

  namespace {
    constexpr double kPI = 3.14159265358979323846;

    //Hard-sphere SANS demo: amorphous SiO2 (rho=2.2 g/cm3) with solid spheres
    //of radius 100 Angstrom (vacuum as solvent). Cross section values marked
    //as regression anchors are anchored to the released plugin behaviour.
    const char * hardsphere_sio2_data =
      "NCMAT v6\n"
      "@DENSITY\n"
      "  2.2 g_per_cm3\n"
      "@DYNINFO\n"
      "  element Si\n"
      "  fraction 0.3333333\n"
      "  type freegas\n"
      "@DYNINFO\n"
      "  element O\n"
      "  fraction 0.6666667\n"
      "  type freegas\n"
      "@CUSTOM_SASCSNS\n"
      "  HardSphere\n"
      "  radius 100\n";

    //DirectLoad demo: trivial I(Q) table over five Q values. In the low-energy
    //limit (neutron wavelength much larger than 1/Q_min) the total SANS cross
    //section must converge to 4*pi*I(Q_min):
    const char * directload_data =
      "NCMAT v6\n"
      "@DENSITY\n"
      "  2.2 g_per_cm3\n"
      "@DYNINFO\n"
      "  element Si\n"
      "  fraction 0.3333333\n"
      "  type freegas\n"
      "@DYNINFO\n"
      "  element O\n"
      "  fraction 0.6666667\n"
      "  type freegas\n"
      "@CUSTOM_SASCSNS\n"
      "  DirectLoad\n"
      "  Q 0.001 0.002 0.004 0.008 0.016\n"
      "  I 1.0 0.9 0.5 0.1 0.01\n";
  }

  void customPluginTest()
  {
    NCPLUGIN_MSG( "Testing plugin" );

    //Register in-memory test data (virtual files, no cleanup needed):
    NC::registerInMemoryFileData( "sascsns_test_hardsphere.ncmat",
                                  std::string( hardsphere_sio2_data ) );
    NC::registerInMemoryFileData( "sascsns_test_directload.ncmat",
                                  std::string( directload_data ) );

    if ( true ) {
      //Check the HardSphere mode (100 Angstrom radius SiO2 spheres):
      auto scat = NC::createScatter( "virtual::sascsns_test_hardsphere.ncmat" );

      //SANS xs must grow monotonically towards low energies:
      double xs_1mev = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-2 } ).dbl();
      double xs_9mev = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-3 } ).dbl();
      double xs_01mev = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-4 } ).dbl();
      NCPLUGIN_MSG( "HardSphere SiO2 R=100AA xs: " << xs_1mev << " barn (E=10meV), "
                    << xs_9mev << " barn (E=1meV), " << xs_01mev << " barn (E=0.1meV)" );
      nc_assert_always( xs_01mev > xs_9mev && xs_9mev > xs_1mev );

      //Regression anchor:
      nc_assert_always( NC::floateq( xs_9mev, 224.2271, 1.0e-4 ) );

      //Scattering must be elastic, with |mu| <= 1:
      auto outcome = scat.sampleScatterIsotropic( NC::NeutronEnergy{ 1e-3 } );
      nc_assert_always( NC::floateq( outcome.ekin.dbl(), 1e-3 ) );
      nc_assert_always( std::abs( outcome.mu.dbl() ) <= 1.0 );

      //The "sans" cfg parameter must disable the plugin (falling back to the
      //standard models, with a much smaller cross section):
      auto scat_off = NC::createScatter( "virtual::sascsns_test_hardsphere.ncmat;sans=0" );
      double xs_off = scat_off.crossSectionIsotropic( NC::NeutronEnergy{ 1e-3 } ).dbl();
      NCPLUGIN_MSG( "sans=0 xs at E=1meV: " << xs_off << " barn (stdlib only)" );
      nc_assert_always( xs_off < 0.1 * xs_9mev );
    }

    if ( true ) {
      //Check the DirectLoad mode. In the limit E -> 0, all events sample Q
      //near Q_min, and the total SANS cross section must converge to
      //4*pi*I(Q_min):
      auto scat = NC::createScatter( "virtual::sascsns_test_directload.ncmat" );
      double xs_lowE = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-15 } ).dbl();
      NCPLUGIN_MSG( "DirectLoad xs at E=1e-15eV: " << xs_lowE << " barn (expected "
                    << 4*kPI << " = 4*pi*I(0))" );
      nc_assert_always( NC::floateq( xs_lowE, 4*kPI, 1.0e-5 ) );

      //Regression anchor:
      double xs_1mev = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-3 } ).dbl();
      nc_assert_always( NC::floateq( xs_1mev, 0.017117, 1.0e-3 ) );
    }

    if ( true ) {
      //Check that the data file shipped in data/ loads and matches the
      //HardSphere in-memory twin:
      auto scat = NC::createScatter( "plugins::SasCSNS/sascsns_sio2_spheres.ncmat" );
      double xs = scat.crossSectionIsotropic( NC::NeutronEnergy{ 1e-3 } ).dbl();
      nc_assert_always( NC::floateq( xs, 224.2271, 1.0e-4 ) );
    }

    NCPLUGIN_MSG( "All tests of plugin were successful!" );
  }
}
