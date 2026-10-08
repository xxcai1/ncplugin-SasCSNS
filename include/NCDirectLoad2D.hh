#ifndef NCPlugin_DirectLoad2D_hh
#define NCPlugin_DirectLoad2D_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)
#include <vector>

namespace NCPluginNamespace {

  ///////////////////////////////////////////////////////////////////////////////
  // DirectLoad2D: tabulated anisotropic SANS (spec: doc/anisotropic_directload2d.tex)
  //
  // The model consumes a 2D intensity I(Qx,Qy) [barn/(atom sr)] on a uniform
  // ascending grid, defined as the in-plane (x,y) components of the elastic
  // scattering vector Q = kf - ki in the material frame. Cross sections and
  // sampling follow section 4.2 / 7.3 of the spec:
  //
  //   sigma(ki,E) = 2 * double-integral over the elastic disc D of
  //                 I / (k^2 sqrt(k^2 - rho'^2)) dA          (two branches)
  //
  // evaluated through the smooth alpha-form (rho = k sin(alpha)):
  //   sigma = int_0^pi int_0^2pi I(disc point) sin(alpha) dphi dalpha
  //
  // Sampling uses per-energy CDFs over the outcome angles (theta,psi) built
  // directly from the image (exact, rejection-free, for beams within
  // kOnAxisTilt of +z), a dilated-proposal + pointwise-mask path for beams
  // tilted by up to kMaxTilt, and a uniform-sphere rejection fallback beyond.
  // All 24 validation checks of script/validate_directload2d.py target these
  // exact formulas.
  ///////////////////////////////////////////////////////////////////////////////

  class DirectLoad2D final : public NC::MoveOnly {
  public:

    //Parse the DirectLoad2D section (raises BadInput on any violation).
    static DirectLoad2D createFromInfo( const NC::Info& info );

    static bool isDirectLoad2D( const NC::Info& info );

    //True if the section is present and requests DirectLoad2D.

    double crossSection( double ekin /*eV*/, const NC::NeutronDirection& dir ) const;

    struct Outcome { double ekin_final, ux, uy, uz; };
    Outcome sample( NC::RNG& rng, double ekin /*eV*/, const NC::NeutronDirection& dir ) const;

    static constexpr double kOnAxisTilt = 0.005;   // rad, exact-CDF tier
    static constexpr double kMaxTilt    = 0.060;   // rad, dilated-proposal tier

  private:

    //Bilinear evaluation, 0 outside the grid:
    double evalI( double qx, double qy ) const;
    double evalIdil( double qx, double qy ) const;

    //Per-energy cumulative sampling table over outcome angles (theta,psi).
    //Two bin layouts (field layout):
    //  theta-uniform (dilated proposals): bin ith covers
    //    theta in [ ith*pi/nth, (ith+1)*pi/nth ], sampled uniformly in theta.
    //  q-uniform (plain CDFs, adaptive resolution): bin ith covers
    //    Q_perp = k*sin(theta) in [ ith*k/nth, (ith+1)*k/nth ], i.e.
    //    theta = asin(Q_perp/k) -- resolution follows the table pitch,
    //    concentrating bins in the forward cone where SANS intensity (and
    //    anisotropy) lives, so narrow patterns are not aliased by fixed
    //    16 mrad bins (the "angular-CDF staircase", doc rung 7); sampled
    //    uniformly in Q_perp within a bin.
    struct AngularCDF {
      double k = 0.0;              // wavenumber [1/Aa]
      double total = 0.0;          // integral of I sin(th) dth dps over sphere
      unsigned nth = 0;            // theta-bin count of this CDF
      bool quniform = false;       // true: Q_perp-uniform bins (see above)
      std::vector<float> cum;      // normalised cumulative cell masses
    };

    void init();//build energy grid, dilated image, CDFs, envelopes, sigma grid
    void buildEnergyData( unsigned ie, double k );//worker: one energy node

    //One-shot runtime warnings (mutable: warned from const methods; plain
    //bools are fine -- NCrystal process objects are single-threaded):
    mutable bool m_warned_sclamp = false;
    mutable bool m_warned_ewindow = false;

    void buildAngularCDF( AngularCDF& cdf, double k, bool dilated ) const;

    //Grid data (I values stored row-major, qy outer loop):
    std::vector<double> m_qx, m_qy;
    std::size_t m_nx = 0, m_ny = 0;
    std::vector<float> m_vals, m_vals_dil;
    double m_qx0 = 0.0, m_qy0 = 0.0, m_invdx = 0.0, m_invdy = 0.0;
    double m_i_max = 0.0;

    //Energy grid (ascending, eV) + per-energy CDFs (plain and dilated):
    std::vector<double> m_ekin;
    std::vector<AngularCDF> m_cdfs, m_cdfs_dil;
    //kNth: theta bins of the dilated proposals and the floor resolution of
    //the plain CDFs. Plain CDFs adapt their theta-bin count to the table
    //pitch: n = 2*k/h (h = finest grid pitch), clamped to
    //[kNth,kNthMaxPlain]; the cap bounds the memory of the 48 plain CDFs
    //(kNthMaxPlain*kNpsi*4B*kNEnergy ~ 38 MB worst case).
    static constexpr unsigned kNth = 192, kNthMaxPlain = 512;
    static constexpr unsigned kNpsi = 384, kNEnergy = 48;

    //Sigma grid over (s,psi_beam,E) with s = k*sin(beta) the in-plane beam
    //shift, multi-linear in all three (s beyond the grid clamps; low energies
    //thereby cover ANY beam tilt, since s <= k there):
    std::vector<double> m_sigma;
    static constexpr unsigned kNS = 41, kNPsi = 41;
    static constexpr double kSMax = 6.948 * kMaxTilt;
    double sigmaLookup( double s, double psi_beam, double ekin_eV ) const;

    //Tier-1 tilt correction envelopes (one per energy node):
    std::vector<double> m_tilt_envelope;

    std::size_t energyNode( double ekin_eV ) const;

    static constexpr double kOfEkinEMax() { return 6.948; } // k of the 100 meV grid top

    //Energy segment (and linear-in-lnE interpolation weight t in [0,1]):
    void energySegment( double ekin_eV, std::size_t& ie0, double& t ) const;
  };

  //The NCrystal process wrapping the model (anisotropic material type):
  class DirectLoad2DScatter final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:
    const char * name() const noexcept override { return NCPLUGIN_NAME_CSTR "DirectLoad2D"; }
    explicit DirectLoad2DScatter( DirectLoad2D&& m ) : m_model(std::move(m)) {}

    NC::CrossSect crossSection( NC::CachePtr&, NC::NeutronEnergy ekin,
                                const NC::NeutronDirection& dir ) const override
    {
      return NC::CrossSect{ m_model.crossSection( ekin.dbl(), dir ) };
    }

    NC::ScatterOutcome sampleScatter( NC::CachePtr&, NC::RNG& rng,
                                      NC::NeutronEnergy ekin,
                                      const NC::NeutronDirection& dir ) const override
    {
      auto o = m_model.sample( rng, ekin.dbl(), dir );
      return { NC::NeutronEnergy{o.ekin_final}, NC::NeutronDirection{o.ux, o.uy, o.uz} };
    }

  private:
    DirectLoad2D m_model;
  };

}

#endif
