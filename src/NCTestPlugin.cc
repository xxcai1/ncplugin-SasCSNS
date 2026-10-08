
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

    //DirectLoad2D demo tables. Constant 3x3 table: in the low-energy limit
    //(disc well inside the table) sigma must be 4*pi*I0 for any beam
    //direction. Peaked table: regression anchor for the full numerics
    //(sigma grid, CDF totals, lnE interpolation).
    const char * directload2d_const_data =
      "NCMAT v7\n"
      "@DENSITY\n"
      "  2.2 g_per_cm3\n"
      "@DYNINFO\n"
      "  element Si\n"
      "  fraction 1.0\n"
      "  type freegas\n"
      "@CUSTOM_SASCSNS\n"
      "  DirectLoad2D\n"
      "  Qx -0.5 0.0 0.5\n"
      "  Qy -0.5 0.0 0.5\n"
      "  I 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n";

    const char * directload2d_peaked_data =
      "NCMAT v7\n"
      "@DENSITY\n"
      "  2.2 g_per_cm3\n"
      "@DYNINFO\n"
      "  element Si\n"
      "  fraction 1.0\n"
      "  type freegas\n"
      "@CUSTOM_SASCSNS\n"
      "  DirectLoad2D\n"
      "  Qx -0.5 0.0 0.5\n"
      "  Qy -0.5 0.0 0.5\n"
      "  I 0.1 0.2 0.1\n"
      "    0.2 0.5 0.2\n"
      "    0.1 0.2 0.1\n";

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
    NC::registerInMemoryFileData( "sascsns_test_directload2d_const.ncmat",
                                  std::string( directload2d_const_data ) );
    NC::registerInMemoryFileData( "sascsns_test_directload2d_peaked.ncmat",
                                  std::string( directload2d_peaked_data ) );

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

    if ( true ) {
      //Check the DirectLoad2D mode (tabulated anisotropic 2D SANS). The
      //constant 3x3 table must give sigma = 4*pi*I0 = 2*pi at low energy,
      //for on-axis, tilted (tier 1), and perpendicular (tier 2 / s-grid)
      //beam directions alike:
      auto scat = NC::createScatter( "virtual::sascsns_test_directload2d_const.ncmat" );
      const double sigma_const = 4*kPI*0.5;
      double xs_z    = scat.crossSection( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{0,0,1} ).dbl();
      double xs_tilt = scat.crossSection( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{0.03,0,0.99955} ).dbl();
      double xs_x    = scat.crossSection( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{1,0,0} ).dbl();
      NCPLUGIN_MSG( "DirectLoad2D const-table xs (E=0.1meV): " << xs_z
                    << " (on axis), " << xs_tilt << " (30mrad), "
                    << xs_x << " (perp); expected " << sigma_const );
      nc_assert_always( NC::floateq( xs_z, sigma_const, 1.0e-4 ) );
      //off-axis rows additionally carry the alpha-quadrature discretisation
      //of the s-grid (midpoint error ~1e-4 on sin(alpha)):
      nc_assert_always( NC::floateq( xs_tilt, sigma_const, 2.0e-4 ) );
      nc_assert_always( NC::floateq( xs_x,    sigma_const, 2.0e-4 ) );

      //Scattering is elastic with unit-modulus outcome directions:
      auto out = scat.sampleScatter( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{0,0,1} );
      nc_assert_always( out.ekin.dbl() == 1e-4 );
      double norm = std::sqrt( out.direction[0]*out.direction[0]
                               + out.direction[1]*out.direction[1]
                               + out.direction[2]*out.direction[2] );
      nc_assert_always( NC::floateq( norm, 1.0, 1.0e-12 ) );

      //The "sans" parameter must disable the plugin here as well:
      auto scat_off = NC::createScatter( "virtual::sascsns_test_directload2d_const.ncmat;sans=0" );
      double xs_off = scat_off.crossSection( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{0,0,1} ).dbl();
      NCPLUGIN_MSG( "DirectLoad2D sans=0 xs at E=0.1meV: " << xs_off << " (stdlib freegas only)" );
      //the stdlib freegas value (7.11 barn, ~1/E rise) must differ clearly
      //from the SANS value (2*pi = 6.283), proving the plugin was bypassed:
      nc_assert_always( !NC::floateq( xs_off, sigma_const, 0.01 ) );

      //Peaked table: regression anchors for the sigma grid / CDF totals /
      //lnE interpolation numerics (values anchored to the released plugin):
      auto scat2 = NC::createScatter( "virtual::sascsns_test_directload2d_peaked.ncmat" );
      double xs_pk   = scat2.crossSection( NC::NeutronEnergy{1e-4}, NC::NeutronDirection{0,0,1} ).dbl();
      double xs_pk2  = scat2.crossSection( NC::NeutronEnergy{5e-4}, NC::NeutronDirection{0,0,1} ).dbl();
      NCPLUGIN_MSG( "DirectLoad2D peaked-table xs: " << xs_pk
                    << " (E=0.1meV), " << xs_pk2 << " (E=0.5meV)" );
      nc_assert_always( NC::floateq( xs_pk,  4.72985, 1.0e-4 ) );
      nc_assert_always( NC::floateq( xs_pk2, 3.09409, 1.0e-4 ) );

      //Peak must scatter forward at k >> table extent (E=5meV gives
      //k=1.55 Angstrom^-1 vs table half-width 0.5: |Q|<=k confines
      //sin(theta) <= 0.5/1.55 = 0.32). Outcomes may land on either elastic
      //branch, hence the absolute value:
      auto out2 = scat2.sampleScatter( NC::NeutronEnergy{5e-3}, NC::NeutronDirection{0,0,1} );
      nc_assert_always( std::abs( out2.direction[2] ) > 0.9 );
    }

    NCPLUGIN_MSG( "All tests of plugin were successful!" );
  }
}
