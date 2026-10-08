#include "NCDirectLoad2D.hh"

#include <atomic>
#include <mutex>
#include <thread>

//Include various utilities from NCrystal's internal header files:
#include "NCrystal/internal/utils/NCString.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

namespace NCPluginNamespace {

  namespace {

    //Constant of the model (spec section 2): E[meV] = eOverK2 * k^2 [1/Aa^2]
    constexpr double eOverK2 = 2.07214;

    inline double kOfEkinE(double ekin_eV)
    {
      return std::sqrt( ( ekin_eV * 1000.0 ) / eOverK2 );
    }

    double str2dblOrThrow( const std::string& s, const char* what )
    {
      double v(0.);
      if ( ! NC::safe_str2dbl( s, v ) )
        NCRYSTAL_THROW2( BadInput, "Invalid numeric value in " << what << ": \"" << s << "\"" );
      return v;
    }

    //1D running maximum with half-width w along one axis of a row-major
    //flattened image (used for the dilated proposal; a square dilation
    //dominates the disc dilation, which is all the tilted-beam mask needs).
    void maxFilterAxis( const std::vector<float>& src, std::vector<float>& dst,
                        std::size_t outer, std::size_t inner, std::size_t w )
    {
      dst.assign( src.size(), 0.0f );
      for ( std::size_t o = 0; o < outer; ++o ) {
        const float* line = &src[ o * inner ];
        float* out = &dst[ o * inner ];
        for ( std::size_t i = 0; i < inner; ++i ) {
          const std::size_t lo = ( i > w ? i - w : 0 );
          const std::size_t hi = std::min( i + w + 1, inner );
          float m = line[ lo ];
          for ( std::size_t j = lo + 1; j < hi; ++j )
            m = std::max( m, line[j] );
          out[i] = m;
        }
      }
    }

  } // namespace

  bool DirectLoad2D::isDirectLoad2D( const NC::Info& info )
  {
    if ( info.countCustomSections( pluginNameUpperCase() ) == 0 )
      return false;
    const auto data = info.getCustomSection( pluginNameUpperCase(), 0 );
    for ( const auto& line : data )
      if ( ! line.empty() && line.at(0) == "DirectLoad2D" )
        return true;
    return false;
  }

  DirectLoad2D DirectLoad2D::createFromInfo( const NC::Info& info )
  {
    const auto data = info.getCustomSection( pluginNameUpperCase(), 0 );

    //Locate the grammar elements: the Qx and Qy axis lines and the I values
    //(free line wrapping after the I keyword, spec section 5.2).
    auto itQx = data.end();
    auto itQy = data.end();
    auto itI  = data.end();
    for ( auto it = data.begin(); it != data.end(); ++it ) {
      if ( it->empty() )
        continue;
      if ( it->at(0) == "Qx" && itQx == data.end() )
        itQx = it;
      else if ( it->at(0) == "Qy" && itQy == data.end() )
        itQy = it;
      else if ( it->at(0) == "I" && itI == data.end() )
        itI = it;
    }
    if ( itQx == data.end() || itQy == data.end() || itI == data.end() )
      NCRYSTAL_THROW( BadInput, "DirectLoad2D section must contain Qx, Qy and I lines" );

    DirectLoad2D m;

    //Axes: at least 2 strictly ascending, uniformly spaced values.
    auto parseAxis = [ &data ]( const decltype(itQx)& it,
                                const char* name,
                                std::vector<double>& ax ) {
      if ( it == data.end() || it->size() < 3 )
        NCRYSTAL_THROW2( BadInput, "DirectLoad2D " << name
                         << " axis needs at least 2 values" );
      ax.reserve( it->size() - 1 );
      for ( std::size_t i = 1; i < it->size(); ++i )
        ax.push_back( str2dblOrThrow( it->at(i), name ) );
      const double d0 = ax[1] - ax[0];
      if ( ! ( d0 > 0.0 ) )
        NCRYSTAL_THROW2( BadInput, "DirectLoad2D " << name
                         << " axis must be strictly ascending" );
      for ( std::size_t i = 2; i < ax.size(); ++i ) {
        if ( ! ( ax[i] > ax[i-1] ) )
          NCRYSTAL_THROW2( BadInput, "DirectLoad2D " << name
                           << " axis must be strictly ascending" );
        if ( std::abs( ( ax[i] - ax[i-1] ) - d0 ) > 1e-6 * std::abs( d0 ) )
          NCRYSTAL_THROW2( BadInput, "DirectLoad2D " << name
                           << " axis must be uniformly spaced" );
      }
    };
    parseAxis( itQx, "Qx", m.m_qx );
    parseAxis( itQy, "Qy", m.m_qy );
    m.m_nx = m.m_qx.size();
    m.m_ny = m.m_qy.size();
    m.m_qx0 = m.m_qx.front();
    m.m_qy0 = m.m_qy.front();
    m.m_invdx = 1.0 / ( m.m_qx.at(1) - m.m_qx.at(0) );
    m.m_invdy = 1.0 / ( m.m_qy.at(1) - m.m_qy.at(0) );

    //I values, row major (qy outer loop), free line wrapping: consume all
    //lines after the I keyword until the next known keyword line.
    m.m_vals.reserve( m.m_nx * m.m_ny );
    //values may start on the I keyword line itself (spec grammar):
    for ( std::size_t t = 1; t < itI->size(); ++t )
      m.m_vals.push_back( static_cast<float>( str2dblOrThrow( itI->at(t), "I values" ) ) );
    for ( auto it = std::next(itI); it != data.end(); ++it ) {
      if ( ! it->empty()
           && ( it->at(0) == "Qx" || it->at(0) == "Qy" || it->at(0) == "I"
                || it->at(0) == "DirectLoad2D" ) )
        break;
      for ( const auto& tok : *it )
        m.m_vals.push_back( static_cast<float>( str2dblOrThrow( tok, "I values" ) ) );
    }
    if ( m.m_vals.size() != m.m_nx * m.m_ny )
      NCRYSTAL_THROW2( BadInput, "DirectLoad2D I section must contain exactly "
                       << ( m.m_nx * m.m_ny ) << " values (" << m.m_nx << "x"
                       << m.m_ny << "), got " << m.m_vals.size() );
    for ( float v : m.m_vals )
      if ( ! ( v >= 0.0f ) )
        NCRYSTAL_THROW( BadInput, "DirectLoad2D I values must be non-negative" );
    if ( *std::max_element( m.m_vals.begin(), m.m_vals.end() ) <= 0.0f )
      NCRYSTAL_THROW( BadInput, "DirectLoad2D I values are all zero" );

    m.init();
    return m;
  }

  void DirectLoad2D::buildEnergyData( unsigned ie, double k )
  {
    //Per-energy CDFs, plain and dilated (spec 7.3):
    buildAngularCDF( m_cdfs[ie], k, false );
    buildAngularCDF( m_cdfs_dil[ie], k, true );

    //Tier-1 tilt envelope: max ratio I(shifted beam)/I(z-hat beam) over the
    //sampling grid, a few tilt azimuths, times a small safety margin:
    {
      const unsigned nth = 96, npsi = 192;
      const double dth = NC::kPi / nth, dps = 2.0 * NC::kPi / npsi;
      const double b1 = kOnAxisTilt;
      double m = 1.0;
      for ( unsigned ia = 0; ia < 8; ++ia ) {
        const double psb = ia * NC::kPi / 4.0;
        const double ax = std::sin( b1 ) * std::cos( psb );
        const double ay = std::sin( b1 ) * std::sin( psb );
        for ( unsigned ith = 0; ith < nth; ++ith ) {
          const double th = ( ith + 0.5 ) * dth;
          for ( unsigned ipsi = 0; ipsi < npsi; ++ipsi ) {
            const double psi = ( ipsi + 0.5 ) * dps;
            const double kfx = std::sin( th ) * std::cos( psi );
            const double kfy = std::sin( th ) * std::sin( psi );
            const double iz = evalI( k * kfx, k * kfy );
            if ( iz <= 0.0 )
              continue;
            const double it = evalI( k * ( kfx - ax ), k * ( kfy - ay ) );
            const double r = it / iz;
            if ( r > m )
              m = r;
          }
        }
      }
      m_tilt_envelope[ie] = 1.01 * m;
    }

    //Sigma grid slice at this energy: (s = k*sin(beta) in [0,kSMax]) x
    //psi_beam in [0,2pi], midpoint quadrature of the smooth alpha-form
    //(spec eq. 11). Multi-linear at run time in (s,psi,lnE); s beyond the
    //grid clamps, which within the design range never happens
    //(s_max = k_top*beta_max covers any beam tilt at E <= E_top).
    {
      constexpr unsigned nalpha = 64, nphi = 64;
      const double dalpha = NC::kPi / nalpha, dphi = 2.0 * NC::kPi / nphi;
      const double ds = kSMax / ( kNS - 1 );
      const double dpsi = 2.0 * NC::kPi / ( kNPsi - 1 );
      double alpha_weight[ nalpha ];
      for ( unsigned ia = 0; ia < nalpha; ++ia )
        alpha_weight[ia] = std::sin( ( ia + 0.5 ) * dalpha ) * dalpha;
      for ( unsigned is_ = 0; is_ < kNS; ++is_ ) {
        const double shift = is_ * ds;
        const double cx0 = -shift;
        for ( unsigned ip = 0; ip < kNPsi; ++ip ) {
          const double psb = ip * dpsi;
          const double cx = cx0 * std::cos( psb ), cy = cx0 * std::sin( psb );
          double acc = 0.0;
          for ( unsigned ial = 0; ial < nalpha; ++ial ) {
            const double al = ( ial + 0.5 ) * dalpha;
            const double ra = k * std::sin( al ), w = alpha_weight[ial];
            for ( unsigned iph = 0; iph < nphi; ++iph ) {
              const double ph = ( iph + 0.5 ) * dphi;
              acc += evalI( cx + ra * std::cos( ph ), cy + ra * std::sin( ph ) ) * w;
            }
          }
          m_sigma[ ( is_ * kNPsi + ip ) * kNEnergy + ie ] = acc * dphi;
        }
      }
    }
  }

  void DirectLoad2D::init()
  {
    m_i_max = *std::max_element( m_vals.begin(), m_vals.end() );

    //Energy grid: 24 log-spaced points in [0.1, 100] meV (spec 7.3). The
    //neighbouring-node spacing keeps sigma- and sampling-grid errors at the
    //per-mille level for smooth tables.
    m_ekin.resize( kNEnergy );
    {
      const double emin = 0.1, emax = 100.0; //meV
      for ( unsigned i = 0; i < kNEnergy; ++i )
        m_ekin[i] = emin * std::pow( emax / emin, double(i) / double( kNEnergy - 1 ) ) * 0.001; //eV
    }

    //Dilated image: square max-filter with half-width = maximal in-plane beam
    //shift k_max*sin(kMaxTilt), in grid cells.
    {
      const double shift = kOfEkinE( m_ekin.back() ) * kMaxTilt;
      const std::size_t wx = std::size_t( std::ceil( shift * m_invdx ) );
      const std::size_t wy = std::size_t( std::ceil( shift * m_invdy ) );
      std::vector<float> tmp;
      maxFilterAxis( m_vals, tmp, m_ny, m_nx, std::min( wx, m_nx - 1 ) );
      maxFilterAxis( tmp, m_vals_dil, m_nx, m_ny, std::min( wy, m_ny - 1 ) );
    }

    //Pre-size everything; the per-energy work below is independent per
    //energy node and runs on worker threads (spec 7.3).
    m_cdfs.resize( kNEnergy );
    m_cdfs_dil.resize( kNEnergy );
    m_tilt_envelope.assign( kNEnergy, 1.0 );
    m_sigma.assign( kNS * kNPsi * kNEnergy, 0.0 );
    {
      std::vector<double> kk( kNEnergy );
      for ( unsigned i = 0; i < kNEnergy; ++i )
        kk[i] = kOfEkinE( m_ekin[i] );
      const unsigned nthr
        = std::min<unsigned>( kNEnergy, std::thread::hardware_concurrency() );
      std::atomic<unsigned> next{ 0 };
      std::exception_ptr eptr;
      std::mutex eptr_mutex;
      auto worker = [&]() {
        for (;;) {
          const unsigned ie = next.fetch_add( 1 );
          if ( ie >= kNEnergy || eptr )
            return;
          try {
            buildEnergyData( ie, kk[ie] );
          } catch ( ... ) {
            std::lock_guard<std::mutex> guard( eptr_mutex );
            if ( !eptr )
              eptr = std::current_exception();
            return;
          }
        }
      };
      std::vector<std::thread> threads;
      threads.reserve( nthr );
      for ( unsigned i = 0; i < nthr; ++i )
        threads.emplace_back( worker );
      for ( auto& t : threads )
        t.join();
      if ( eptr )
        std::rethrow_exception( eptr );
    }

    //Overwrite the s = 0 row with the CDF totals: the on-axis sigma is then
    //the very same discretised integral the samplers normalise to, making
    //xs/sampler consistency exact at the energy nodes (and preserved by the
    //lnE interpolation, since the mixture sampler uses the same weights).
    //This also sidesteps quadrature accuracy loss on table-truncated discs.
    for ( unsigned ie = 0; ie < kNEnergy; ++ie )
      for ( unsigned ip = 0; ip < kNPsi; ++ip )
        m_sigma[ ( 0 * kNPsi + ip ) * kNEnergy + ie ] = m_cdfs[ie].total;
  }

  void DirectLoad2D::buildAngularCDF( AngularCDF& cdf, double k, bool dilated ) const
  {
    cdf.k = k;
    cdf.cum.assign( kNth * kNpsi, 0.0f );
    const double dth = NC::kPi / kNth, dps = 2.0 * NC::kPi / kNpsi;
    double tot = 0.0;
    for ( unsigned ith = 0; ith < kNth; ++ith ) {
      const double th = ( ith + 0.5 ) * dth;
      const double sth = std::sin( th );
      const double w = sth * dth * dps;
      for ( unsigned ipsi = 0; ipsi < kNpsi; ++ipsi ) {
        const double psi = ( ipsi + 0.5 ) * dps;
        const double qx = k * sth * std::cos( psi );
        const double qy = k * sth * std::sin( psi );
        const double v = dilated ? evalIdil( qx, qy ) : evalI( qx, qy );
        tot += v * w;
        cdf.cum[ ith * kNpsi + ipsi ] = static_cast<float>( tot );
      }
    }
    cdf.total = tot;
    if ( tot > 0.0 ) {
      const double inv = 1.0 / tot;
      for ( auto& x : cdf.cum )
        x = static_cast<float>( x * inv );
      cdf.cum.back() = 1.0f;
    }
  }

  double DirectLoad2D::evalI( double qx, double qy ) const
  {
    const double fx = ( qx - m_qx0 ) * m_invdx;
    const double fy = ( qy - m_qy0 ) * m_invdy;
    if ( fx < 0.0 || fy < 0.0 || fx > double(m_nx - 1) || fy > double(m_ny - 1) )
      return 0.0;
    const std::size_t ix = std::min( std::size_t(fx), m_nx - 2 );
    const std::size_t iy = std::min( std::size_t(fy), m_ny - 2 );
    const double tx = fx - ix, ty = fy - iy;
    const float* row0 = &m_vals[ iy * m_nx + ix ];
    const float* row1 = row0 + m_nx;
    return ( row0[0] * ( 1.0 - tx ) + row0[1] * tx ) * ( 1.0 - ty )
         + ( row1[0] * ( 1.0 - tx ) + row1[1] * tx ) * ty;
  }

  double DirectLoad2D::evalIdil( double qx, double qy ) const
  {
    const double fx = ( qx - m_qx0 ) * m_invdx;
    const double fy = ( qy - m_qy0 ) * m_invdy;
    if ( fx < 0.0 || fy < 0.0 || fx > double(m_nx - 1) || fy > double(m_ny - 1) )
      return 0.0;
    const std::size_t ix = std::min( std::size_t(fx), m_nx - 2 );
    const std::size_t iy = std::min( std::size_t(fy), m_ny - 2 );
    const double tx = fx - ix, ty = fy - iy;
    const float* row0 = &m_vals_dil[ iy * m_nx + ix ];
    const float* row1 = row0 + m_nx;
    return ( row0[0] * ( 1.0 - tx ) + row0[1] * tx ) * ( 1.0 - ty )
         + ( row1[0] * ( 1.0 - tx ) + row1[1] * tx ) * ty;
  }

  std::size_t DirectLoad2D::energyNode( double ekin_eV ) const
  {
    auto it = std::lower_bound( m_ekin.begin(), m_ekin.end(), ekin_eV );
    if ( it == m_ekin.end() )
      return m_ekin.size() - 1;
    if ( it == m_ekin.begin() )
      return 0;
    const std::size_t i = it - m_ekin.begin();
    return ( ekin_eV - m_ekin[i-1] < m_ekin[i] - ekin_eV ) ? i - 1 : i;
  }

  void DirectLoad2D::energySegment( double ekin_eV, std::size_t& ie0, double& t ) const
  {
    //grid is log-uniform in energy: locate the segment and the linear-in-lnE
    //weight; clamp outside the grid range.
    if ( ekin_eV <= m_ekin.front() ) {
      ie0 = 0; t = 0.0; return;
    }
    if ( ekin_eV >= m_ekin.back() ) {
      ie0 = m_ekin.size() - 2; t = 1.0; return;
    }
    auto it = std::lower_bound( m_ekin.begin(), m_ekin.end(), ekin_eV );
    const std::size_t i = it - m_ekin.begin();
    ie0 = i - 1;
    t = ( std::log( ekin_eV ) - std::log( m_ekin[ie0] ) )
      / ( std::log( m_ekin[i] ) - std::log( m_ekin[ie0] ) );
  }

  double DirectLoad2D::sigmaLookup( double s, double psi_beam, double ekin_eV ) const
  {
    while ( psi_beam < 0.0 )
      psi_beam += 2.0 * NC::kPi;
    while ( psi_beam >= 2.0 * NC::kPi )
      psi_beam -= 2.0 * NC::kPi;
    const double dpsi = 2.0 * NC::kPi / ( kNPsi - 1 );
    const double ds = kSMax / ( kNS - 1 );
    const double fs = std::min( s, kSMax ) / ds;
    const double fp = psi_beam / dpsi;
    const std::size_t is_ = std::min( std::size_t(fs), std::size_t(kNS - 2) );
    const std::size_t ip = std::min( std::size_t(fp), std::size_t(kNPsi - 2) );
    const double ts = std::min( fs - is_, 1.0 );
    const double tp = std::min( fp - ip, 1.0 );
    std::size_t ie; double te;
    energySegment( ekin_eV, ie, te );
    auto at = [&]( std::size_t sidx, std::size_t pidx, std::size_t eidx ) -> const double&
    {
      return m_sigma[ ( sidx * kNPsi + pidx ) * kNEnergy + eidx ];
    };
    double sum = 0.0;
    const double w[2] = { 1.0 - te, te };
    for ( unsigned de = 0; de < 2; ++de ) {
      const double s00 = at( is_, ip, ie + de ), s10 = at( is_ + 1, ip, ie + de );
      const double s01 = at( is_, ip + 1, ie + de ), s11 = at( is_ + 1, ip + 1, ie + de );
      sum += w[de] * ( ( s00 * ( 1 - ts ) + s10 * ts ) * ( 1 - tp )
                       + ( s01 * ( 1 - ts ) + s11 * ts ) * tp );
    }
    return sum;
  }

  double DirectLoad2D::crossSection( double ekin, const NC::NeutronDirection& dir ) const
  {
    const double beta = std::acos( std::clamp( dir[2], -1.0, 1.0 ) );
    const double psi_beam = std::atan2( dir[1], dir[0] );
    return sigmaLookup( kOfEkinE( ekin ) * std::sin( beta ), psi_beam, ekin );
  }

  DirectLoad2D::Outcome DirectLoad2D::sample( NC::RNG& rng, double ekin,
                                              const NC::NeutronDirection& dir ) const
  {
    Outcome o;
    o.ekin_final = ekin;//elastic
    o.ux = dir[0]; o.uy = dir[1]; o.uz = dir[2];

    const double k = kOfEkinE( ekin );
    std::size_t ie0; double te;
    energySegment( ekin, ie0, te );
    const double ax = dir[0], ay = dir[1];
    const double beta = std::acos( std::clamp( dir[2], -1.0, 1.0 ) );
    //Node mixture weights, consistent with the lnE-linear sigma interpolation:
    const std::size_t nodes[2] = { ie0, ie0 + 1 };
    const double wmix[2] = { 1.0 - te, te };

    //Sample a direction from a per-energy angular CDF:
    auto sampleAngles = [ &rng ]( const AngularCDF& cdf,
                                  double& kfx, double& kfy, double& kfz ) {
      const double u = rng.generate();
      auto it = std::upper_bound( cdf.cum.begin(), cdf.cum.end(), u );
      std::size_t idx = ( it == cdf.cum.end() ? cdf.cum.size() - 1 : it - cdf.cum.begin() );
      const double dth = NC::kPi / kNth, dps = 2.0 * NC::kPi / kNpsi;
      const std::size_t ith = idx / kNpsi, ipsi = idx % kNpsi;
      const double th = ( ith + rng.generate() ) * dth;
      const double psi = ( ipsi + rng.generate() ) * dps;
      const double st = std::sin( th );
      kfx = st * std::cos( psi );
      kfy = st * std::sin( psi );
      kfz = std::cos( th );
    };

    if ( beta <= kOnAxisTilt ) {
      //Tier 1: exact on-axis CDFs of the two adjacent energy nodes, mixed
      //with the sigma interpolation weights (+ ratio mask with scanned
      //envelope for small residual tilts). No rejection at beta == 0.
      if ( ! ( m_cdfs[ nodes[0] ].total > 0.0 || m_cdfs[ nodes[1] ].total > 0.0 ) )
        return o;//nothing to scatter from at this energy
      const double M0 = m_tilt_envelope[ nodes[0] ], M1 = m_tilt_envelope[ nodes[1] ];
      for ( int itry = 0; itry < 10000; ++itry ) {
        const bool second = rng.generate() >= wmix[0];
        const AngularCDF& cdf = m_cdfs[ nodes[ second ? 1 : 0 ] ];
        const double M = second ? M1 : M0;
        if ( ! ( cdf.total > 0.0 ) )
          continue;
        double kfx, kfy, kfz;
        sampleAngles( cdf, kfx, kfy, kfz );
        if ( beta < 1.0e-12 ) {
          o.ux = kfx; o.uy = kfy; o.uz = kfz;
          return o;
        }
        const double km = cdf.k;
        const double iz = evalI( km * kfx, km * kfy );
        if ( ! ( iz > 0.0 ) )
          continue;//outside proposal support
        const double it = evalI( km * ( kfx - ax ), km * ( kfy - ay ) );
        if ( rng.generate() * M * iz <= it ) {
          o.ux = kfx; o.uy = kfy; o.uz = kfz;
          return o;
        }
      }
    }
    else if ( beta <= kMaxTilt ) {
      //Tier 2: dilated proposal + pointwise mask (bounded by 1 by the
      //dilation argument, spec 7.3).
      if ( ! ( m_cdfs_dil[ nodes[0] ].total > 0.0 || m_cdfs_dil[ nodes[1] ].total > 0.0 ) )
        return o;
      for ( int itry = 0; itry < 10000; ++itry ) {
        const bool second = rng.generate() >= wmix[0];
        const AngularCDF& cdf = m_cdfs_dil[ nodes[ second ? 1 : 0 ] ];
        if ( ! ( cdf.total > 0.0 ) )
          continue;
        double kfx, kfy, kfz;
        sampleAngles( cdf, kfx, kfy, kfz );
        const double km = cdf.k;
        const double idil = evalIdil( km * kfx, km * kfy );
        if ( ! ( idil > 0.0 ) )
          continue;
        const double it = evalI( km * ( kfx - ax ), km * ( kfy - ay ) );
        if ( rng.generate() * idil <= it ) {
          o.ux = kfx; o.uy = kfy; o.uz = kfz;
          return o;
        }
      }
    }

    //Fallback (beyond design tilt, or tiers exhausted): uniform-sphere
    //rejection on the plain image. Slow but exact.
    if ( ! ( m_i_max > 0.0 ) )
      return o;
    for ( int itry = 0; itry < 100000; ++itry ) {
      const double z = 2.0 * rng.generate() - 1.0;
      const double phi = 2.0 * NC::kPi * rng.generate();
      const double r = std::sqrt( std::max( 0.0, 1.0 - z * z ) );
      const double kfx = r * std::cos( phi ), kfy = r * std::sin( phi ), kfz = z;
      const double it = evalI( k * ( kfx - ax ), k * ( kfy - ay ) );
      if ( rng.generate() * m_i_max <= it ) {
        o.ux = kfx; o.uy = kfy; o.uz = kfz;
        return o;
      }
    }
    return o;//pathological; leave unchanged rather than loop forever
  }

} // namespace NCPluginNamespace
